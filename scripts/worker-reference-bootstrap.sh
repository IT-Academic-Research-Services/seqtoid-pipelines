#!/bin/bash
set -uo pipefail

LOG=/var/log/userdata-sandbox.log
exec > >(tee -a "$LOG") 2>&1

echo "=== AL2023 user-data started at $(date -Is) ==="
log() { echo "[$(date -Is)] $*"; }

# ============================================================
# Worker identity / registration metadata
# ============================================================
# These describe what this instance is intended to be.
# They are registration metadata, not AWS EC2 tags.
WORKER_ROLE="${WORKER_ROLE:-seqtoid-nr-cpu-worker}"
WORKER_BACKEND="${WORKER_BACKEND:-cpu}"
WORKER_REFERENCE_SET="${WORKER_REFERENCE_SET:-phase2}"
WORKER_ENVIRONMENT="${WORKER_ENVIRONMENT:-dev}"

WORKER_STATUS_DIR="/efs/workers"

# ---------- helpers ----------
dnf_retry() {
  local tries=10
  local i
  for i in $(seq 1 "$tries"); do
    if dnf -y install "$@"; then
      return 0
    fi
    log "dnf failed (try $i/$tries) - waiting for RPM lock"
    sleep 15
  done
  log "ERROR: dnf still failing after $tries tries: $*"
  return 1
}

# ---------- EC2 metadata ----------
get_imds_token() {
  curl -s -X PUT \
    "http://169.254.169.254/latest/api/token" \
    -H "X-aws-ec2-metadata-token-ttl-seconds: 21600" \
    || true
}

imds_get() {
  local token="$1"
  local path="$2"

  curl -s \
    -H "X-aws-ec2-metadata-token: $token" \
    "http://169.254.169.254/latest/meta-data/${path}" \
    || true
}

# ---------- THP ----------
setup_thp() {
  log "THP"

  if [ -w /sys/kernel/mm/transparent_hugepage/enabled ]; then
    echo madvise > /sys/kernel/mm/transparent_hugepage/enabled
    echo defer+madvise > /sys/kernel/mm/transparent_hugepage/defrag 2>/dev/null || true
    echo never > /sys/kernel/mm/transparent_hugepage/shmem_enabled 2>/dev/null || true
    echo 1 > /sys/kernel/mm/transparent_hugepage/khugepaged/defrag 2>/dev/null || true
  fi
}

# ---------- packages ----------
install_packages() {
  log "Packages - critical"

  if ! dnf_retry mdadm xfsprogs nfs-utils amazon-ssm-agent libatomic; then
    log "ERROR: critical package installation failed"
    return 1
  fi

  systemctl enable --now amazon-ssm-agent || true

  log "Packages - secondary"
  dnf_retry nvme-cli pigz git tar wget python3-pip || true
}

# ---------- s5cmd ----------
# Pin the worker image to a known s5cmd release. s5cmd is the reference
# staging tool for large S3 -> NVMe transfers; do not fall back to aws s3 cp.
install_s5cmd() {
  local version="2.3.0"
  local url="https://github.com/peak/s5cmd/releases/download/v${version}/s5cmd_${version}_Linux-64bit.tar.gz"
  local tmpdir="/tmp/s5cmd-install"
  local archive="${tmpdir}/s5cmd.tar.gz"

  log "Installing s5cmd ${version}"

  rm -rf "$tmpdir"
  mkdir -p "$tmpdir"

  if ! curl -fL --retry 5 --retry-delay 2 --connect-timeout 15 \
      -o "$archive" "$url"; then
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

  if ! /usr/local/bin/s5cmd version; then
    log "ERROR: s5cmd installation verification failed"
    return 1
  fi

  log "s5cmd installed: $(/usr/local/bin/s5cmd version | head -n 1)"
  return 0
}

# ---------- /scratch helper (runs on every boot via systemd) ----------
# The RAID is built ONLY from EC2 instance-store NVMe devices.
# There is no EBS RAID setup here.
#
# NEVER put /dev/md0 in fstab. Instance-store dies on stop/start;
# fstab would risk an emergency boot.
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

# Let NVMe nodes appear
sleep 8

# 1) Try assemble anything with existing superblocks (plain reboot)
mdadm --assemble --scan 2>/dev/null || true

# 2) Prefer an existing unmounted RAID0
EXISTING=$(lsblk -nr -o NAME,TYPE,MOUNTPOINT \
  | awk '$2=="raid0" && $3=="" {print "/dev/"$1; exit}')

if [ -n "${EXISTING:-}" ]; then
  RAID="$EXISTING"
  echo "Using existing RAID $RAID"
else
  # Discover EC2 instance-store NVMe devices.
  # These are local ephemeral NVMe devices, NOT EBS volumes.
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

  # 3) Assemble by member list if superblocks exist (reboot)
  if mdadm --assemble "$RAID" "${DEVS[@]}" 2>/dev/null; then
    echo "Assembled $RAID from members"
  else
    # 4) Stop/start wiped store -> create fresh RAID0
    mdadm --stop "$RAID" 2>/dev/null || true

    echo "Creating RAID0 $RAID (${#DEVS[@]} devices, chunk 256k)"
    mdadm --create --verbose --force "$RAID" \
      --level=0 \
      --chunk=256 \
      --raid-devices="${#DEVS[@]}" \
      "${DEVS[@]}"
  fi
fi

# Wait for array to go clean
sleep 2

if [ ! -b "$RAID" ]; then
  # Fallback: first unmounted raid0
  RAID=$(lsblk -nr -o NAME,TYPE,MOUNTPOINT \
    | awk '$2=="raid0" && $3=="" {print "/dev/"$1; exit}')

  if [ -z "${RAID:-}" ]; then
    echo "ERROR: no RAID block device available"
    cat /proc/mdstat || true
    exit 1
  fi
fi

# Format only if no filesystem (new array after stop/start)
if ! blkid "$RAID" 2>/dev/null | grep -q 'TYPE='; then
  N=$(mdadm --detail "$RAID" 2>/dev/null \
    | awk '/Raid Devices/ {print $4}')

  N=${N:-4}

  echo "mkfs.xfs -f -L scratch -d su=256k,sw=$N -l size=128m $RAID"

  mkfs.xfs -f \
    -L scratch \
    -d su=256k,sw="$N" \
    -l size=128m \
    "$RAID"
fi

mount \
  -o noatime,nodiratime,logbufs=8,logbsize=256k,largeio,inode64,swalloc \
  "$RAID" \
  "$MOUNT"

chmod 1777 "$MOUNT"

echo 8192 > \
  /sys/block/$(basename "$RAID")/queue/read_ahead_kb \
  2>/dev/null || true

echo "/scratch ready: $(df -h "$MOUNT" | tail -1)"
cat /proc/mdstat || true
SCRIPTEOF

  chmod 755 /usr/local/sbin/setup-scratch.sh

  # Remove any stale fstab line that would break boot
  if grep -qsE '[[:space:]]/scratch[[:space:]]' /etc/fstab; then
    log "Removing /scratch from fstab (unsafe for instance-store RAID)"
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

  log "scratch-setup.service enabled"
}

# ---------- EFS ----------
setup_efs() {
  local EFS_ID="${EFS_ID:-fs-040f78f17e586fc54}"
  local AWS_REGION="${AWS_REGION:-us-west-2}"
  local EFS_DNS="${EFS_ID}.efs.${AWS_REGION}.amazonaws.com"
  local MOUNT=/efs

  log "EFS ${EFS_ID} (${EFS_DNS}) -> $MOUNT"

  mkdir -p "$MOUNT"

  if mountpoint -q "$MOUNT"; then
    log "EFS already mounted"
    return 0
  fi

  if mount -t nfs4 \
    -o nfsvers=4.1,rsize=1048576,wsize=1048576,hard,timeo=600,retrans=2,noresvport,_netdev \
    "${EFS_DNS}:/" "$MOUNT"; then

    log "EFS mounted"
  else
    log "ERROR: EFS mount failed"
    return 1
  fi

  if ! grep -qs "$EFS_DNS" /etc/fstab; then
    echo "${EFS_DNS}:/ $MOUNT nfs4 nfsvers=4.1,rsize=1048576,wsize=1048576,hard,timeo=600,retrans=2,noresvport,_netdev 0 0" \
      >> /etc/fstab
  fi

  chmod 1777 "$MOUNT"
}

# ---------- worker registration/status ----------
write_worker_status() {
  local state="$1"
  local worker_ready="$2"
  local reference_status="$3"
  local reason="${4:-}"

  log "Writing worker registration/status: state=${state} worker_ready=${worker_ready} reference=${reference_status}"

  local token
  local instance_id
  local instance_type
  local private_ip
  local placement
  local availability_zone
  local ami_id
  local hostname
  local scratch_size
  local scratch_used
  local scratch_available
  local scratch_device
  local efs_status
  local status_dir
  local status_file
  local tmp_file
  local now

  token=$(get_imds_token)

  instance_id=$(imds_get "$token" "instance-id")
  instance_type=$(imds_get "$token" "instance-type")
  private_ip=$(imds_get "$token" "local-ipv4")
  placement=$(imds_get "$token" "placement/availability-zone")
  availability_zone="$placement"
  ami_id=$(imds_get "$token" "ami-id")

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
    scratch_size=0
    scratch_used=0
    scratch_available=0
    scratch_device="not-mounted"
  fi

  if mountpoint -q /efs; then
    efs_status="mounted"
  else
    efs_status="not-mounted"
  fi

  status_dir="${WORKER_STATUS_DIR}/${instance_id}"
  status_file="${status_dir}/status.json"
  tmp_file="${status_file}.tmp"

  mkdir -p "$status_dir"

  now=$(date -Is)

  python3 - "$tmp_file" \
    "$now" \
    "$state" \
    "$worker_ready" \
    "$reason" \
    "$instance_id" \
    "$instance_type" \
    "$private_ip" \
    "$availability_zone" \
    "$ami_id" \
    "$hostname" \
    "$WORKER_ROLE" \
    "$WORKER_BACKEND" \
    "$WORKER_REFERENCE_SET" \
    "$WORKER_ENVIRONMENT" \
    "$scratch_size" \
    "$scratch_used" \
    "$scratch_available" \
    "$scratch_device" \
    "$efs_status" \
    "$reference_status" << 'PYEOF'
import json
import sys

(
    output_path,
    timestamp,
    state,
    worker_ready,
    reason,
    instance_id,
    instance_type,
    private_ip,
    availability_zone,
    ami_id,
    hostname,
    worker_role,
    backend,
    reference_set,
    environment,
    scratch_size,
    scratch_used,
    scratch_available,
    scratch_device,
    efs_status,
    reference_status,
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
        "efs": {
            "status": efs_status,
        },
    },

    "reference": {
        "status": reference_status,
        "reference_set": reference_set,
    },

    "registered_at": timestamp,
    "last_update": timestamp,
}

if reason:
    status["reason"] = reason

with open(output_path, "w", encoding="utf-8") as fh:
    json.dump(status, fh, indent=2, sort_keys=True)
    fh.write("\n")
PYEOF

  if [ ! -f "$tmp_file" ]; then
    log "ERROR: failed to create worker status file"
    return 1
  fi

  chmod 644 "$tmp_file"
  mv -f "$tmp_file" "$status_file"

  log "Wrote worker status: $status_file"
  log "Worker state: ${state}"
  log "Worker ready: ${worker_ready}"
  log "Worker identity: ${instance_id} (${instance_type})"
  log "Worker role: ${WORKER_ROLE}"
  log "Worker backend: ${WORKER_BACKEND}"
  log "Reference set: ${WORKER_REFERENCE_SET} (${reference_status})"
  log "Environment: ${WORKER_ENVIRONMENT}"

  if [ -n "$reason" ]; then
    log "Worker reason: ${reason}"
  fi
}

# ---------- readiness transitions ----------
set_worker_state() {
  local state="$1"
  local worker_ready="$2"
  local reference_status="$3"
  local reason="${4:-}"

  if ! write_worker_status \
      "$state" \
      "$worker_ready" \
      "$reference_status" \
      "$reason"; then
    log "WARNING: failed to publish worker state=${state}"
  fi
}

# ---------- main ----------
# A worker is not READY merely because EC2 reports it as running.
# Reference preparation and validation must complete before READY.

setup_thp

if ! install_packages; then
  log "ERROR: critical package installation failed; worker cannot be used"
  set_worker_state \
    "FAILED" \
    "false" \
    "not_checked" \
    "critical package installation failed"
  exit 1
fi

if ! install_scratch_helper; then
  log "ERROR: failed to install scratch setup helper; worker cannot be used"
  set_worker_state \
    "FAILED" \
    "false" \
    "not_checked" \
    "failed to install scratch setup helper"
  exit 1
fi

# Let systemd be the authoritative mechanism for preparing /scratch.
# This reconstructs the filesystem/mount; it does NOT preserve instance-store data.
# EC2 instance-store NVMe contents are ephemeral and may disappear on stop/start.

if ! systemctl start scratch-setup.service; then
  set_worker_state \
    "FAILED" \
    "false" \
    "not_checked" \
    "instance-store /scratch setup service failed"

  log "ERROR: scratch-setup.service failed; worker cannot be used"
  exit 1
fi

# A successful Type=oneshot + RemainAfterExit service should remain active.
if ! systemctl is-active --quiet scratch-setup.service; then
  set_worker_state \
    "FAILED" \
    "false" \
    "not_checked" \
    "scratch-setup.service did not become active"

  log "ERROR: scratch-setup.service is not active after start"
  systemctl status scratch-setup.service --no-pager || true
  exit 1
fi

if ! setup_efs; then
  set_worker_state \
    "FAILED" \
    "false" \
    "not_checked" \
    "EFS setup failed"

  log "ERROR: EFS setup failed; worker cannot be used"
  exit 1
fi

mkdir -p /efs/workers

if ! install_s5cmd; then
  set_worker_state \
    "FAILED" \
    "false" \
    "not_checked" \
    "s5cmd installation failed"

  log "ERROR: s5cmd installation failed; worker cannot be used"
  exit 1
fi

set_worker_state \
  "BOOTING" \
  "false" \
  "not_checked"

set_worker_state \
  "PREPARING_REFERENCE" \
  "false" \
  "preparing"

REFERENCE_PREPARER="/usr/local/bin/seqtoid_prepare_reference.py"

if [ ! -f "$REFERENCE_PREPARER" ]; then
  set_worker_state \
    "FAILED" \
    "false" \
    "failed" \
    "reference preparation utility not found: $REFERENCE_PREPARER"

  log "ERROR: reference preparation utility not found: $REFERENCE_PREPARER"
  exit 1
fi

case "$WORKER_BACKEND" in
  cpu)
    REFERENCE_BACKEND="mmseqs-cpu"
    ;;

  mmseqs-cpu)
    REFERENCE_BACKEND="mmseqs-cpu"
    ;;

  mmseqs-gpu)
    REFERENCE_BACKEND="mmseqs-gpu"
    ;;

  diamond)
    REFERENCE_BACKEND="diamond"
    ;;

  *)
    set_worker_state \
      "FAILED" \
      "false" \
      "failed" \
      "unsupported worker backend=$WORKER_BACKEND"

    log "ERROR: unsupported worker backend=$WORKER_BACKEND"
    exit 1
    ;;
esac

log "Preparing reference for backend=$REFERENCE_BACKEND"

REFERENCE_VERSION=$(python3 "$REFERENCE_PREPARER" \
  --backend "$REFERENCE_BACKEND" \
  --reference-set "$WORKER_REFERENCE_SET") || {
    set_worker_state \
      "FAILED" \
      "false" \
      "failed" \
      "reference preparation failed for backend=$REFERENCE_BACKEND"

    log "ERROR: reference preparation failed for backend=$REFERENCE_BACKEND"
    exit 1
  }

case "$REFERENCE_BACKEND" in
  mmseqs-cpu)
    REFERENCE_DIR="/scratch/refs/mmseqs"
    REFERENCE_DB="$REFERENCE_DIR/nrcleanDB"
    ;;

  mmseqs-gpu)
    REFERENCE_DIR="/scratch/refs/mmseqs-gpu"
    REFERENCE_DB="$REFERENCE_DIR/nrcleanDB_gpu"
    ;;

  diamond)
    REFERENCE_DIR="/scratch/refs/diamond"
    REFERENCE_DB="$REFERENCE_DIR/diamond_07_22_2026.dmnd"
    ;;
esac

if [ -z "$REFERENCE_VERSION" ] || \
   [ ! -f "$REFERENCE_DIR/.reference_version" ] || \
   [ ! -f "$REFERENCE_DIR/.reference_manifest.json" ] || \
   [ ! -f "$REFERENCE_DB" ]; then

  set_worker_state \
    "FAILED" \
    "false" \
    "failed" \
    "reference validation outputs missing"

  log "ERROR: reference validation outputs missing under $REFERENCE_DIR"
  exit 1
fi

RECORDED_VERSION=$(tr -d '[:space:]' < "$REFERENCE_DIR/.reference_version")

if [ "$REFERENCE_VERSION" != "$RECORDED_VERSION" ]; then
  set_worker_state \
    "FAILED" \
    "false" \
    "failed" \
    "reference version verification failed"

  log "ERROR: reference version mismatch: preparer=$REFERENCE_VERSION local=$RECORDED_VERSION"
  exit 1
fi

log "Reference prepared and validated:"
log "  backend=$REFERENCE_BACKEND"
log "  reference_set=$WORKER_REFERENCE_SET"
log "  version=$REFERENCE_VERSION"
log "  directory=$REFERENCE_DIR"
log "  database=$REFERENCE_DB"

set_worker_state \
  "READY" \
  "true" \
  "validated"

mkdir -p /home/ec2-user/{programs,venv} /dev/shm/workspace

chown -R ec2-user:ec2-user \
  /home/ec2-user/programs \
  /home/ec2-user/venv \
  || true

log "=== user-data finished at $(date -Is) ==="

log "Verify:"
log "  df -h /scratch /efs"
log "  systemctl is-enabled scratch-setup.service"
log "  cat /proc/mdstat"

log "Worker status files:"
find /efs/workers \
  -maxdepth 2 \
  -type f \
  -name 'status.json' \
  -print
