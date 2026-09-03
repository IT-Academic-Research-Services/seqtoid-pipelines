#!/bin/bash
set -uo pipefail

LOG=/var/log/userdata-gpu.log
exec > >(tee -a "$LOG") 2>&1

echo "=== AL2023 MMseqs GPU worker bootstrap started at $(date -Is) ==="
log() { echo "[$(date -Is)] $*"; }

# GPU worker identity / registration metadata
WORKER_ROLE="${WORKER_ROLE:-seqtoid-nr-mmseqs-gpu-worker}"
WORKER_BACKEND="${WORKER_BACKEND:-}"
WORKER_REFERENCE_SET="${WORKER_REFERENCE_SET:-phase2}"
WORKER_ENVIRONMENT="${WORKER_ENVIRONMENT:-dev}"
WORKER_STATUS_DIR="/efs/workers"

# Helpers

dnf_retry() {
  local tries=10
  local i
  for i in $(seq 1 "$tries"); do
    if dnf -y install "$@"; then
      return 0
    fi
    log "dnf failed (try $i/$tries) - waiting"
    sleep 15
  done
  log "ERROR: dnf still failing after $tries tries: $*"
  return 1
}

get_imds_token() {
  curl -sS -X PUT \
    "http://169.254.169.254/latest/api/token" \
    -H "X-aws-ec2-metadata-token-ttl-seconds: 21600" || true
}

imds_get() {
  local token="$1"
  local path="$2"
  curl -sS -H "X-aws-ec2-metadata-token: $token" \
    "http://169.254.169.254/latest/meta-data/${path}" || true
}

validate_worker_configuration() {
  if [ -z "$WORKER_BACKEND" ]; then
    log "ERROR: WORKER_BACKEND was not supplied by the GPU launch wrapper"
    return 1
  fi
  if [ "$WORKER_BACKEND" != "mmseqs-gpu" ]; then
    log "ERROR: GPU bootstrap requires WORKER_BACKEND=mmseqs-gpu, got $WORKER_BACKEND"
    return 1
  fi
  return 0
}

setup_thp() {
  log "THP"
  if [ -w /sys/kernel/mm/transparent_hugepage/enabled ]; then
    echo madvise > /sys/kernel/mm/transparent_hugepage/enabled
    echo defer+madvise > /sys/kernel/mm/transparent_hugepage/defrag 2>/dev/null || true
    echo never > /sys/kernel/mm/transparent_hugepage/shmem_enabled 2>/dev/null || true
    echo 1 > /sys/kernel/mm/transparent_hugepage/khugepaged/defrag 2>/dev/null || true
  fi
}

install_packages() {
  log "Packages - critical"
  if ! dnf_retry mdadm xfsprogs nfs-utils amazon-ssm-agent libatomic python3; then
    log "ERROR: critical package installation failed"
    return 1
  fi
  systemctl enable --now amazon-ssm-agent || true

  log "Packages - secondary"
  dnf_retry nvme-cli pigz git tar wget python3-pip || true
}

install_s5cmd() {
  local version="2.3.0"
  local url="https://github.com/peak/s5cmd/releases/download/v${version}/s5cmd_${version}_Linux-64bit.tar.gz"
  local tmpdir="/tmp/s5cmd-install"
  local archive="${tmpdir}/s5cmd.tar.gz"

  log "Installing s5cmd ${version}"
  rm -rf "$tmpdir"
  mkdir -p "$tmpdir"

  if ! curl -fL --retry 5 --retry-delay 2 --connect-timeout 15 -o "$archive" "$url"; then
    log "ERROR: failed to download s5cmd ${version}"
    rm -rf "$tmpdir"
    return 1
  fi
  if ! tar -xzf "$archive" -C "$tmpdir" s5cmd; then
    log "ERROR: failed to extract s5cmd ${version}"
    rm -rf "$tmpdir"
    return 1
  fi
  install -m 0755 "$tmpdir/s5cmd" /usr/local/bin/s5cmd
  rm -rf "$tmpdir"
  /usr/local/bin/s5cmd version >/dev/null || {
    log "ERROR: s5cmd verification failed"
    return 1
  }
  log "s5cmd installed: $(/usr/local/bin/s5cmd version | head -n 1)"
}

validate_gpu_runtime() {
  log "Validating NVIDIA GPU runtime"
  if ! command -v nvidia-smi >/dev/null 2>&1; then
    log "ERROR: nvidia-smi not found; install the NVIDIA compute driver before running this bootstrap"
    return 1
  fi
  nvidia-smi

  if ! nvidia-smi -L >/dev/null 2>&1; then
    log "ERROR: NVIDIA driver cannot enumerate GPUs"
    return 1
  fi

  if ! command -v mmseqs >/dev/null 2>&1; then
    log "ERROR: mmseqs executable not found"
    return 1
  fi

  log "MMseqs version:"
  mmseqs version || true

  if ! mmseqs search --help 2>&1 | grep -q -- "--gpu"; then
    log "ERROR: installed mmseqs does not advertise --gpu support"
    return 1
  fi

  log "NVIDIA + MMseqs GPU runtime validated"
}

install_scratch_helper() {
  log "Installing /usr/local/sbin/setup-scratch.sh + systemd unit"

  cat > /usr/local/sbin/setup-scratch.sh << 'SCRIPTEOF'
#!/bin/bash
set -uo pipefail
MOUNT=/scratch
RAID=/dev/md0
LOG=/var/log/setup-scratch.log
exec >> "$LOG" 2>&1

echo "=== setup-scratch $(date -Is) ==="
mkdir -p "$MOUNT"
if mountpoint -q "$MOUNT"; then
  echo "/scratch already mounted"
  df -h "$MOUNT" | tail -1
  exit 0
fi
sleep 8
mdadm --assemble --scan 2>/dev/null || true
EXISTING=$(lsblk -nr -o NAME,TYPE,MOUNTPOINT | awk '$2=="raid0" && $3=="" {print "/dev/"$1; exit}')
if [ -n "${EXISTING:-}" ]; then
  RAID="$EXISTING"
  echo "Using existing RAID $RAID"
else
  mapfile -t DEVS < <(
    lsblk -dn -o NAME,MODEL \
      | awk 'tolower($0) ~ /instance storage/ {print "/dev/"$1}' \
      | sort
  )
  if [ ${#DEVS[@]} -eq 0 ]; then
    echo "ERROR: no instance-store NVMe found"
    lsblk -d -o NAME,SIZE,TYPE,MODEL
    exit 1
  fi
  echo "Instance-store devices: ${DEVS[*]}"
  if mdadm --assemble "$RAID" "${DEVS[@]}" 2>/dev/null; then
    echo "Assembled $RAID from members"
  else
    mdadm --stop "$RAID" 2>/dev/null || true
    echo "Creating RAID0 $RAID (${#DEVS[@]} devices, chunk 256k)"
    mdadm --create --verbose --force "$RAID" \
      --level=0 --chunk=256 --raid-devices="${#DEVS[@]}" "${DEVS[@]}"
  fi
fi
sleep 2
if [ ! -b "$RAID" ]; then
  RAID=$(lsblk -nr -o NAME,TYPE,MOUNTPOINT \
    | awk '$2=="raid0" && $3=="" {print "/dev/"$1; exit}')
  if [ -z "${RAID:-}" ]; then
    echo "ERROR: no RAID block device available"
    cat /proc/mdstat || true
    exit 1
  fi
fi
if ! blkid "$RAID" 2>/dev/null | grep -q 'TYPE='; then
  N=$(mdadm --detail "$RAID" 2>/dev/null | awk '/Raid Devices/ {print $4}')
  N=${N:-4}
  mkfs.xfs -f -L scratch -d su=256k,sw="$N" -l size=128m "$RAID"
fi
mount -o noatime,nodiratime,logbufs=8,logbsize=256k,largeio,inode64,swalloc "$RAID" "$MOUNT"
chmod 1777 "$MOUNT"
echo 8192 > /sys/block/$(basename "$RAID")/queue/read_ahead_kb 2>/dev/null || true
echo "/scratch ready: $(df -h "$MOUNT" | tail -1)"
cat /proc/mdstat || true
SCRIPTEOF

  chmod 755 /usr/local/sbin/setup-scratch.sh

  if grep -qsE '[[:space:]]/scratch[[:space:]]' /etc/fstab; then
    log "Removing /scratch from fstab"
    sed -i '\#[[:space:]]/scratch[[:space:]]#d' /etc/fstab
  fi

  cat > /etc/systemd/system/scratch-setup.service << 'EOF'
[Unit]
Description=Assemble/create instance-store RAID0 and mount /scratch
DefaultDependencies=no
Wants=systemd-udevd.service
After=systemd-udevd.service local-fs-pre.target
Before=local-fs.target
Conflicts=shutdown.target

[Service]
Type=oneshot
RemainAfterExit=yes
ExecStart=/usr/local/sbin/setup-scratch.sh
TimeoutStartSec=300
StandardOutput=journal+console
StandardError=journal+console

[Install]
WantedBy=local-fs.target multi-user.target
EOF

  systemctl daemon-reload
  systemctl enable scratch-setup.service
}

setup_efs() {
  local EFS_ID="${EFS_ID:-fs-040f78f17e586fc54}"
  local AWS_REGION="${AWS_REGION:-us-west-2}"
  local EFS_DNS="${EFS_ID}.efs.${AWS_REGION}.amazonaws.com"
  local MOUNT=/efs
  log "EFS ${EFS_ID} (${EFS_DNS}) -> ${MOUNT}"
  mkdir -p "$MOUNT"
  if mountpoint -q "$MOUNT"; then
    log "EFS already mounted"
    return 0
  fi
  if mount -t nfs4 -o nfsvers=4.1,rsize=1048576,wsize=1048576,hard,timeo=600,retrans=2,noresvport,_netdev "${EFS_DNS}:/" "$MOUNT"; then
    log "EFS mounted"
  else
    log "ERROR: EFS mount failed"
    return 1
  fi
  if ! grep -qs "$EFS_DNS" /etc/fstab; then
    echo "${EFS_DNS}:/ $MOUNT nfs4 nfsvers=4.1,rsize=1048576,wsize=1048576,hard,timeo=600,retrans=2,noresvport,_netdev 0 0" >> /etc/fstab
  fi
  chmod 1777 "$MOUNT"
}

write_worker_status() {
  local state="$1"
  local worker_ready="$2"
  local reference_status="$3"
  local reason="${4:-}"
  local token instance_id instance_type private_ip availability_zone ami_id hostname
  local scratch_size scratch_used scratch_available scratch_device efs_status status_dir status_file tmp_file now

  token=$(get_imds_token)
  instance_id=$(imds_get "$token" instance-id)
  instance_type=$(imds_get "$token" instance-type)
  private_ip=$(imds_get "$token" local-ipv4)
  availability_zone=$(imds_get "$token" placement/availability-zone)
  ami_id=$(imds_get "$token" ami-id)
  hostname=$(hostname -s)

  instance_id="${instance_id:-unknown}"
  instance_type="${instance_type:-unknown}"
  private_ip="${private_ip:-unknown}"
  availability_zone="${availability_zone:-unknown}"
  ami_id="${ami_id:-unknown}"

  if mountpoint -q /scratch; then
    scratch_size=$(df -B1 --output=size /scratch | tail -1 | tr -d ' ')
    scratch_used=$(df -B1 --output=used /scratch | tail -1 | tr -d ' ')
    scratch_available=$(df -B1 --output=avail /scratch | tail -1 | tr -d ' ')
    scratch_device=$(findmnt -n -o SOURCE /scratch || echo unknown)
  else
    scratch_size=0; scratch_used=0; scratch_available=0; scratch_device=not-mounted
  fi

  if mountpoint -q /efs; then efs_status=mounted; else efs_status=not-mounted; fi

  status_dir="${WORKER_STATUS_DIR}/${instance_id}"
  status_file="${status_dir}/status.json"
  tmp_file="${status_file}.tmp"
  mkdir -p "$status_dir"
  now=$(date -Is)

  python3 - "$tmp_file" "$now" "$state" "$worker_ready" "$reason" \
    "$instance_id" "$instance_type" "$private_ip" "$availability_zone" "$ami_id" "$hostname" \
    "$WORKER_ROLE" "$WORKER_BACKEND" "$WORKER_REFERENCE_SET" "$WORKER_ENVIRONMENT" \
    "$scratch_size" "$scratch_used" "$scratch_available" "$scratch_device" "$efs_status" "$reference_status" << 'PYEOF'
import json
import sys
(
    output_path, timestamp, state, worker_ready, reason,
    instance_id, instance_type, private_ip, availability_zone, ami_id, hostname,
    worker_role, backend, reference_set, environment,
    scratch_size, scratch_used, scratch_available, scratch_device, efs_status, reference_status,
) = sys.argv[1:]
status = {
    "state": state,
    "worker_ready": worker_ready.lower() == "true",
    "instance": {
        "instance_id": instance_id,
        "instance_type": instance_type,
        "private_ip": private_ip,
        "availability_zone": availability_zone,
        "ami_id": ami_id,
        "hostname": hostname,
    },
    "worker": {
        "role": worker_role,
        "backend": backend,
        "reference_set": reference_set,
        "environment": environment,
    },
    "storage": {
        "scratch": {
            "mounted": scratch_device != "not-mounted",
            "device": scratch_device,
            "size_bytes": int(scratch_size),
            "used_bytes": int(scratch_used),
            "available_bytes": int(scratch_available),
        },
        "efs": {"status": efs_status},
    },
    "reference": {"status": reference_status, "reference_set": reference_set},
    "registered_at": timestamp,
    "last_update": timestamp,
}
if reason:
    status["reason"] = reason
with open(output_path, "w", encoding="utf-8") as fh:
    json.dump(status, fh, indent=2, sort_keys=True)
    fh.write("\n")
PYEOF

  chmod 644 "$tmp_file"
  mv -f "$tmp_file" "$status_file"
  log "Wrote worker status: $status_file"
  log "Worker state: ${state}"
  log "Worker ready: ${worker_ready}"
  log "Worker backend: ${WORKER_BACKEND}"
  log "Reference: ${reference_status} (${WORKER_REFERENCE_SET})"
  [ -n "$reason" ] && log "Worker reason: ${reason}"
}

set_worker_state() {
  if ! write_worker_status "$1" "$2" "$3" "${4:-}"; then
    log "WARNING: failed to publish worker state=$1"
  fi
}

# Main
setup_thp
validate_worker_configuration || exit 1
install_packages || exit 1
install_scratch_helper || exit 1
systemctl start scratch-setup.service || { log "ERROR: scratch setup failed"; exit 1; }
systemctl is-active --quiet scratch-setup.service || { log "ERROR: scratch setup inactive"; exit 1; }
setup_efs || exit 1
mkdir -p /efs/workers

set_worker_state BOOTING false not_checked

install_s5cmd || {
  set_worker_state FAILED false not_checked "s5cmd installation failed"
  exit 1
}

set_worker_state PREPARING_REFERENCE false preparing

validate_gpu_runtime || {
  set_worker_state FAILED false failed "NVIDIA/MMseqs GPU runtime validation failed"
  exit 1
}

REFERENCE_PREPARER="/usr/local/bin/seqtoid_prepare_reference_gpu.py"
REFERENCE_DIR="/scratch/refs/mmseqs-gpu"
REFERENCE_DB="${REFERENCE_DIR}/nrcleanDB_gpu"

if [ ! -f "$REFERENCE_PREPARER" ]; then
  set_worker_state FAILED false failed "reference preparation utility not found: $REFERENCE_PREPARER"
  log "ERROR: reference preparation utility not found: $REFERENCE_PREPARER"
  exit 1
fi

REFERENCE_VERSION=$(python3 "$REFERENCE_PREPARER" --backend mmseqs-gpu --reference-set "$WORKER_REFERENCE_SET") || {
  set_worker_state FAILED false failed "GPU reference preparation failed"
  exit 1
}

if [ -z "$REFERENCE_VERSION" ] || \
   [ ! -f "${REFERENCE_DIR}/.reference_version" ] || \
   [ ! -f "${REFERENCE_DIR}/.reference_manifest.json" ] || \
   [ ! -f "$REFERENCE_DB" ]; then
  set_worker_state FAILED false failed "GPU reference validation outputs missing"
  log "ERROR: reference validation outputs missing under $REFERENCE_DIR"
  exit 1
fi

RECORDED_VERSION=$(tr -d '[:space:]' < "${REFERENCE_DIR}/.reference_version")
if [ "$REFERENCE_VERSION" != "$RECORDED_VERSION" ]; then
  set_worker_state FAILED false failed "GPU reference version verification failed"
  exit 1
fi

log "GPU reference prepared and validated: version=${REFERENCE_VERSION} db=${REFERENCE_DB}"
set_worker_state READY true validated

mkdir -p /home/ec2-user/{programs,venv} /dev/shm/workspace
chown -R ec2-user:ec2-user /home/ec2-user/programs /home/ec2-user/venv || true

log "=== MMseqs GPU worker bootstrap finished at $(date -Is) ==="
log "Verify: df -h /scratch /efs"
log "Verify: nvidia-smi"
log "Verify: mmseqs version"
log "Verify: cat /efs/workers/*/status.json"
