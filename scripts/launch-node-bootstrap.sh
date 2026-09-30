#!/bin/bash
set -uo pipefail

LOG=/var/log/seqtoid-launch-bootstrap.log
exec > >(tee -a "$LOG") 2>&1

echo "=== SeqToID AL2023 launch-node bootstrap started at $(date -Is) ==="

log() {
  echo "[$(date -Is)] $*"
}

# ============================================================
# Launch-node configuration
# ============================================================
#
# This bootstrap is for the SeqToID pipeline launch node.
#
# It is NOT a distributed NR worker bootstrap.
#
# Responsibilities:
#   - configure THP
#   - ensure required OS/runtime packages
#   - ensure SSM is available
#   - install/verify pinned s5cmd
#   - install/enable the instance-store /scratch helper
#   - assemble/create and mount /scratch
#   - mount EFS at /efs
#   - establish:
#
#         /scratch/refs
#         /refs -> /scratch/refs
#
#   - invoke the launch-node reference preparation utility
#
# The reference preparation utility is responsible for synchronizing
# the canonical launch-node reference set from S3.
#
# Launch-node references include everything below:
#
#   s3://seqtoid-public-references/phase2/refs/
#
# EXCEPT exactly:
#
#   nr_clean/cpu/*
#   nr_clean/gpu/*
#   diamond/*
#
# In particular, the launch node DOES require:
#
#   nr_clean/nr_clean.fa
#   nr_clean/nr_clean_index.fst
#
# because those are used by the non-distributed portions of the
# short-read pipeline.
#
# /scratch is EC2 instance-store storage and is ephemeral.
# A stop/start may therefore require the complete reference cache
# to be downloaded again.
# ============================================================

REFERENCE_SET="${REFERENCE_SET:-phase2}"

REFERENCE_ROOT="/scratch/refs"
REFERENCE_LINK="/refs"

REFERENCE_PREPARER="/usr/local/bin/seqtoid_prepare_reference_launch.py"

EFS_ID="${EFS_ID:-fs-040f78f17e586fc54}"
AWS_REGION="${AWS_REGION:-us-west-2}"


# ============================================================
# Helpers
# ============================================================

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


# ============================================================
# Transparent Huge Pages
# ============================================================

setup_thp() {
  log "Configuring Transparent Huge Pages"

  if [ -w /sys/kernel/mm/transparent_hugepage/enabled ]; then
    echo madvise \
      > /sys/kernel/mm/transparent_hugepage/enabled

    echo defer+madvise \
      > /sys/kernel/mm/transparent_hugepage/defrag \
      2>/dev/null || true

    echo never \
      > /sys/kernel/mm/transparent_hugepage/shmem_enabled \
      2>/dev/null || true

    echo 1 \
      > /sys/kernel/mm/transparent_hugepage/khugepaged/defrag \
      2>/dev/null || true
  fi

  log "THP configuration complete"
}


# ============================================================
# Packages
# ============================================================

install_packages() {
  log "Ensuring critical launch-node infrastructure packages"

  if ! dnf_retry \
      mdadm \
      xfsprogs \
      nfs-utils \
      amazon-ssm-agent \
      libatomic; then

    log "ERROR: critical package installation failed"
    return 1
  fi

  systemctl enable --now amazon-ssm-agent || true

  log "Ensuring secondary launch-node infrastructure packages"

  dnf_retry \
    nvme-cli \
    pigz \
    git \
    tar \
    wget \
    python3-pip \
    || true

  return 0
}


# ============================================================
# s5cmd
# ============================================================
#
# Pin the launch-node environment to the same qualified s5cmd
# release used by the worker infrastructure.
#
# s5cmd is used for large S3 -> instance-store transfers.
# ============================================================

install_s5cmd() {
  local version="2.3.0"
  local url
  local tmpdir
  local archive

  url="https://github.com/peak/s5cmd/releases/download/v${version}/s5cmd_${version}_Linux-64bit.tar.gz"

  tmpdir="/tmp/s5cmd-install"
  archive="${tmpdir}/s5cmd.tar.gz"

  # If the correct version is already installed, leave it alone.
  if command -v s5cmd >/dev/null 2>&1; then
    if s5cmd version 2>/dev/null \
        | grep -q "^v${version}"; then

      log "s5cmd ${version} already installed"
      return 0
    fi
  fi

  log "Installing s5cmd ${version}"

  rm -rf "$tmpdir"
  mkdir -p "$tmpdir"

  if ! curl \
      -fL \
      --retry 5 \
      --retry-delay 2 \
      --connect-timeout 15 \
      -o "$archive" \
      "$url"; then

    log "ERROR: failed to download s5cmd ${version}"
    rm -rf "$tmpdir"
    return 1
  fi

  if ! tar \
      -xzf "$archive" \
      -C "$tmpdir" \
      s5cmd; then

    log "ERROR: failed to extract s5cmd ${version}"
    rm -rf "$tmpdir"
    return 1
  fi

  install \
    -m 0755 \
    "$tmpdir/s5cmd" \
    /usr/local/bin/s5cmd

  rm -rf "$tmpdir"

  if ! /usr/local/bin/s5cmd version; then
    log "ERROR: s5cmd installation verification failed"
    return 1
  fi

  log "s5cmd installed: $(/usr/local/bin/s5cmd version | head -n 1)"

  return 0
}


# ============================================================
# /scratch helper
# ============================================================
#
# Runs on every boot via systemd.
#
# The RAID is built ONLY from EC2 instance-store NVMe devices.
#
# There is no EBS RAID setup here.
#
# NEVER put /dev/md0 in /etc/fstab.
#
# Instance-store can disappear on stop/start. A persistent fstab
# entry for /dev/md0 could therefore prevent a future boot.
#
# IMPORTANT:
#
# This launch-node version also creates /scratch/refs with the
# correct ownership after /scratch has been mounted. This prevents
# the root:root /scratch/refs problem observed after rebuilding
# ephemeral scratch.
# ============================================================

install_scratch_helper() {
  log "Installing launch-node /usr/local/sbin/setup-scratch.sh"

  cat > /usr/local/sbin/setup-scratch.sh << 'SCRIPTEOF'
#!/bin/bash
set -uo pipefail

MOUNT="/scratch"
REFERENCE_ROOT="/scratch/refs"
RAID="/dev/md0"

LOG="/var/log/setup-scratch.log"

exec >> "$LOG" 2>&1

echo "=== setup-scratch $(date -Is) ==="

mkdir -p "$MOUNT"


# ------------------------------------------------------------
# Existing mount
# ------------------------------------------------------------

if mountpoint -q "$MOUNT"; then
  echo "/scratch already mounted"

  # The reference directory is part of the launch-node boot
  # contract and must always be writable by ec2-user.
  mkdir -p "$REFERENCE_ROOT"
  chown ec2-user:ec2-user "$REFERENCE_ROOT"
  chmod 755 "$REFERENCE_ROOT"

  df -h "$MOUNT" | tail -1

  exit 0
fi


# ------------------------------------------------------------
# Wait for EC2 NVMe devices
# ------------------------------------------------------------

sleep 8


# ------------------------------------------------------------
# Try to assemble an existing array
# ------------------------------------------------------------
#
# This is the normal reboot case when instance-store contents
# still exist.
# ------------------------------------------------------------

mdadm --assemble --scan 2>/dev/null || true


# ------------------------------------------------------------
# Prefer an existing unmounted RAID0
# ------------------------------------------------------------

EXISTING=$(
  lsblk -nr -o NAME,TYPE,MOUNTPOINT \
    | awk '$2=="raid0" && $3=="" {print "/dev/"$1; exit}'
)

if [ -n "${EXISTING:-}" ]; then

  RAID="$EXISTING"

  echo "Using existing RAID $RAID"

else

  # ----------------------------------------------------------
  # Discover EC2 instance-store NVMe devices
  # ----------------------------------------------------------

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


  # ----------------------------------------------------------
  # Assemble by member list if possible
  # ----------------------------------------------------------

  if mdadm \
      --assemble \
      "$RAID" \
      "${DEVS[@]}" \
      2>/dev/null; then

    echo "Assembled $RAID from members"

  else

    # --------------------------------------------------------
    # Instance-store was lost: create a new RAID0
    # --------------------------------------------------------

    mdadm --stop "$RAID" 2>/dev/null || true

    echo "Creating RAID0 $RAID (${#DEVS[@]} devices, chunk 256k)"

    mdadm \
      --create \
      --verbose \
      --force \
      "$RAID" \
      --level=0 \
      --chunk=256 \
      --raid-devices="${#DEVS[@]}" \
      "${DEVS[@]}"
  fi
fi


# ------------------------------------------------------------
# Wait briefly for array
# ------------------------------------------------------------

sleep 2


# ------------------------------------------------------------
# Resolve RAID block device if necessary
# ------------------------------------------------------------

if [ ! -b "$RAID" ]; then

  RAID=$(
    lsblk -nr -o NAME,TYPE,MOUNTPOINT \
      | awk '$2=="raid0" && $3=="" {print "/dev/"$1; exit}'
  )

  if [ -z "${RAID:-}" ]; then

    echo "ERROR: no RAID block device available"

    cat /proc/mdstat || true

    exit 1
  fi
fi


# ------------------------------------------------------------
# Format only when no filesystem exists
# ------------------------------------------------------------

if ! blkid "$RAID" 2>/dev/null \
    | grep -q 'TYPE='; then

  N=$(
    mdadm --detail "$RAID" 2>/dev/null \
      | awk '/Raid Devices/ {print $4}'
  )

  N="${N:-4}"

  echo "mkfs.xfs -f -L scratch -d su=256k,sw=$N -l size=128m $RAID"

  mkfs.xfs \
    -f \
    -L scratch \
    -d su=256k,sw="$N" \
    -l size=128m \
    "$RAID"
fi


# ------------------------------------------------------------
# Mount
# ------------------------------------------------------------

mount \
  -o noatime,nodiratime,logbufs=8,logbsize=256k,largeio,inode64,swalloc \
  "$RAID" \
  "$MOUNT"


# /scratch itself remains a shared temporary workspace.
chmod 1777 "$MOUNT"


# ------------------------------------------------------------
# Launch-node reference root
# ------------------------------------------------------------

mkdir -p "$REFERENCE_ROOT"

chown \
  ec2-user:ec2-user \
  "$REFERENCE_ROOT"

chmod \
  755 \
  "$REFERENCE_ROOT"


# ------------------------------------------------------------
# RAID read-ahead
# ------------------------------------------------------------

echo 8192 \
  > "/sys/block/$(basename "$RAID")/queue/read_ahead_kb" \
  2>/dev/null || true


echo "/scratch ready: $(df -h "$MOUNT" | tail -1)"

echo "/scratch/refs:"
ls -ld "$REFERENCE_ROOT"

cat /proc/mdstat || true

SCRIPTEOF

  chmod \
    755 \
    /usr/local/sbin/setup-scratch.sh


  # ----------------------------------------------------------
  # Remove unsafe stale /scratch fstab entry
  # ----------------------------------------------------------

  if grep \
      -qsE \
      '[[:space:]]/scratch[[:space:]]' \
      /etc/fstab; then

    log "Removing /scratch from fstab"

    sed \
      -i \
      '\#[[:space:]]/scratch[[:space:]]#d' \
      /etc/fstab
  fi


  # ----------------------------------------------------------
  # systemd service
  # ----------------------------------------------------------

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


# ============================================================
# EFS
# ============================================================

setup_efs() {
  local EFS_DNS
  local MOUNT

  EFS_DNS="${EFS_ID}.efs.${AWS_REGION}.amazonaws.com"
  MOUNT="/efs"

  log "EFS ${EFS_ID} (${EFS_DNS}) -> ${MOUNT}"

  mkdir -p "$MOUNT"

  if mountpoint -q "$MOUNT"; then

    log "EFS already mounted"

  else

    if mount \
        -t nfs4 \
        -o nfsvers=4.1,rsize=1048576,wsize=1048576,hard,timeo=600,retrans=2,noresvport,_netdev \
        "${EFS_DNS}:/" \
        "$MOUNT"; then

      log "EFS mounted"

    else

      log "ERROR: EFS mount failed"

      return 1
    fi
  fi


  if ! grep \
      -qs \
      "$EFS_DNS" \
      /etc/fstab; then

    echo "${EFS_DNS}:/ $MOUNT nfs4 nfsvers=4.1,rsize=1048576,wsize=1048576,hard,timeo=600,retrans=2,noresvport,_netdev 0 0" \
      >> /etc/fstab
  fi

  chmod 1777 "$MOUNT"

  return 0
}


# ============================================================
# /refs launch-node path
# ============================================================

setup_reference_path() {
  log "Configuring launch-node reference path"

  if ! mountpoint -q /scratch; then
    log "ERROR: /scratch is not mounted"
    return 1
  fi

  mkdir -p "$REFERENCE_ROOT"

  chown \
    ec2-user:ec2-user \
    "$REFERENCE_ROOT"

  chmod \
    755 \
    "$REFERENCE_ROOT"


  # ----------------------------------------------------------
  # /refs must either:
  #
  #   not exist
  #
  # or:
  #
  #   already be the correct symlink
  #
  # We deliberately refuse to replace a real /refs directory,
  # because doing so could destroy or hide reference data.
  # ----------------------------------------------------------

  if [ -L "$REFERENCE_LINK" ]; then

    local target

    target=$(
      readlink -f "$REFERENCE_LINK" \
        2>/dev/null || true
    )

    if [ "$target" != "$REFERENCE_ROOT" ]; then

      log "ERROR: $REFERENCE_LINK points to unexpected target: $target"

      return 1
    fi

  elif [ -e "$REFERENCE_LINK" ]; then

    log "ERROR: $REFERENCE_LINK exists and is not a symlink"

    return 1

  else

    ln \
      -s \
      "$REFERENCE_ROOT" \
      "$REFERENCE_LINK"
  fi


  # ----------------------------------------------------------
  # Explicit write test
  # ----------------------------------------------------------

  if ! runuser \
      -u ec2-user \
      -- touch "$REFERENCE_ROOT/.seqtoid-write-test"; then

    log "ERROR: ec2-user cannot write to $REFERENCE_ROOT"

    return 1
  fi

  rm -f "$REFERENCE_ROOT/.seqtoid-write-test"


  log "Launch-node reference path ready:"
  log "  $REFERENCE_LINK -> $(readlink -f "$REFERENCE_LINK")"

  ls -ld \
    "$REFERENCE_LINK" \
    "$REFERENCE_ROOT"

  return 0
}


# ============================================================
# Validate installed launch-node software
# ============================================================
#
# This is deliberately a runtime validation step, not an installer.
#
# The launch-node AMI should already contain the scientific
# software that was compiled/qualified on the target CPU family.
#
# Add/remove entries here if the final pipeline software inventory
# changes before the AMI is frozen.
# ============================================================

validate_launch_software() {
  local required=(
    seqtoid-pipelines
    fastp
    kallisto
    bowtie2
    hisat2
    minimap2
    samtools
    bcftools
    s5cmd
  )

  local executable

  log "Validating launch-node software"

  for executable in "${required[@]}"; do

    if ! command \
        -v "$executable" \
        >/dev/null 2>&1; then

      log "ERROR: required launch-node executable is missing: $executable"

      return 1
    fi

    log "  found $executable: $(command -v "$executable")"
  done


  # Runtime linker checks for native HTSlib consumers.

  if ldd "$(command -v samtools)" \
      | grep -q 'not found'; then

    log "ERROR: samtools has unresolved shared libraries"

    ldd "$(command -v samtools)"

    return 1
  fi

  if ldd "$(command -v bcftools)" \
      | grep -q 'not found'; then

    log "ERROR: bcftools has unresolved shared libraries"

    ldd "$(command -v bcftools)"

    return 1
  fi


  # bcftools must be able to enumerate its plugins.
  # This catches the symbol-visibility failure encountered during
  # the launch-node AMI build.

  if ! bcftools plugin -l \
      >/tmp/seqtoid-bcftools-plugins.txt \
      2>/tmp/seqtoid-bcftools-plugins.err; then

    log "ERROR: bcftools plugin validation failed"

    cat /tmp/seqtoid-bcftools-plugins.err || true

    rm -f \
      /tmp/seqtoid-bcftools-plugins.txt \
      /tmp/seqtoid-bcftools-plugins.err

    return 1
  fi

  if grep \
      -q \
      'undefined symbol' \
      /tmp/seqtoid-bcftools-plugins.err; then

    log "ERROR: bcftools plugins contain unresolved symbols"

    cat /tmp/seqtoid-bcftools-plugins.err

    rm -f \
      /tmp/seqtoid-bcftools-plugins.txt \
      /tmp/seqtoid-bcftools-plugins.err

    return 1
  fi

  rm -f \
    /tmp/seqtoid-bcftools-plugins.txt \
    /tmp/seqtoid-bcftools-plugins.err


  log "Launch-node software validation passed"

  return 0
}


# ============================================================
# Reference preparation
# ============================================================

prepare_references() {
  log "Preparing launch-node references"

  if [ ! -f "$REFERENCE_PREPARER" ]; then

    log "ERROR: launch-node reference preparation utility not found:"
    log "  $REFERENCE_PREPARER"

    return 1
  fi

  if [ ! -x "$REFERENCE_PREPARER" ]; then

    log "ERROR: launch-node reference preparation utility is not executable:"
    log "  $REFERENCE_PREPARER"

    return 1
  fi


  if ! "$REFERENCE_PREPARER" \
      --reference-set "$REFERENCE_SET" \
      --local-root "$REFERENCE_ROOT"; then

    log "ERROR: launch-node reference preparation failed"

    return 1
  fi


  log "Launch-node reference preparation completed"

  return 0
}


# ============================================================
# Main
# ============================================================

setup_thp


if ! install_packages; then

  log "ERROR: critical package installation failed"

  exit 1
fi


if ! install_s5cmd; then

  log "ERROR: s5cmd installation failed"

  exit 1
fi


if ! install_scratch_helper; then

  log "ERROR: failed to install scratch setup helper"

  exit 1
fi


# ------------------------------------------------------------
# systemd is authoritative for /scratch
# ------------------------------------------------------------

if ! systemctl start scratch-setup.service; then

  log "ERROR: scratch-setup.service failed"

  systemctl status \
    scratch-setup.service \
    --no-pager \
    || true

  exit 1
fi


if ! systemctl \
    is-active \
    --quiet \
    scratch-setup.service; then

  log "ERROR: scratch-setup.service is not active after start"

  systemctl status \
    scratch-setup.service \
    --no-pager \
    || true

  exit 1
fi


if ! mountpoint -q /scratch; then

  log "ERROR: /scratch is not mounted after scratch setup"

  exit 1
fi


# ------------------------------------------------------------
# EFS
# ------------------------------------------------------------

if ! setup_efs; then

  log "ERROR: EFS setup failed"

  exit 1
fi


# ------------------------------------------------------------
# /refs
# ------------------------------------------------------------

if ! setup_reference_path; then

  log "ERROR: launch-node reference path setup failed"

  exit 1
fi


# ------------------------------------------------------------
# Scientific software
# ------------------------------------------------------------

if ! validate_launch_software; then

  log "ERROR: launch-node software validation failed"

  exit 1
fi


# ------------------------------------------------------------
# Reference cache
# ------------------------------------------------------------

if ! prepare_references; then

  log "ERROR: launch-node reference preparation failed"

  exit 1
fi


# ============================================================
# Final validation
# ============================================================

if ! mountpoint -q /scratch; then

  log "ERROR: final validation: /scratch is not mounted"

  exit 1
fi


if ! mountpoint -q /efs; then

  log "ERROR: final validation: /efs is not mounted"

  exit 1
fi


if [ "$(readlink -f "$REFERENCE_LINK" 2>/dev/null || true)" \
     != "$REFERENCE_ROOT" ]; then

  log "ERROR: final validation: $REFERENCE_LINK is incorrect"

  exit 1
fi


if ! runuser \
    -u ec2-user \
    -- test \
    -w "$REFERENCE_ROOT"; then

  log "ERROR: final validation: $REFERENCE_ROOT is not writable by ec2-user"

  exit 1
fi


log "============================================================"
log "SeqToID launch node READY"
log "============================================================"

log "Reference set: $REFERENCE_SET"
log "Reference root: $REFERENCE_ROOT"
log "/refs target: $(readlink -f "$REFERENCE_LINK")"
log "Scratch: $(df -h /scratch | tail -1)"
log "EFS: $(findmnt -n -o SOURCE /efs || echo unknown)"

log "=== SeqToID AL2023 launch-node bootstrap finished at $(date -Is) ==="

exit 0
