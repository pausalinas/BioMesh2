#!/usr/bin/env bash
set -o pipefail

SCRIPT_NAME=$(basename "$0")
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEFAULT_CONFIG="${SCRIPT_DIR}/config/default.conf"
VERSION="1.0.0"

print_banner() {
    cat <<'EOF'
============================================================
   ██████  ██ ██████  ███    ███ ███████ ███████ ██   ██
   ██   ██ ██ ██   ██ ████  ████ ██      ██      ██   ██
   ██   ██ ██ ██████  ██ ████ ██ █████   ███████ ███████
   ██   ██ ██ ██   ██ ██  ██  ██ ██           ██ ██   ██
   ██████  ██ ██   ██ ██      ██ ███████ ███████ ██   ██
                     BioMesh2 CLI Runner
============================================================
EOF
}

usage() {
    cat <<EOF
Usage: ${SCRIPT_NAME} [options] <input.pdb>

Options:
  -i, --input FILE           Input PDB file (positional also accepted)
  -o, --output FILE          Output mesh file (default: mesh.vtk)
  -v, --voxel-size VALUE     Voxel size in Angstroms (default: 1.0)
  -p, --padding VALUE        Bounding box padding in Angstroms (default: 0.0)
  -f, --format FORMAT        Output format (vtk|gid|json|txt) (default: vtk)
      --filter TYPE          Filter type (none|all|protein-only|no-water|custom)
  -c, --config FILE          INI-style config file
      --generate-config      Print default config and exit
      --batch LIST           Comma-separated list of PDB files to process
  -V, --verbose              Verbose output
  -q, --quiet                Minimal output
      --version              Show version
  -h, --help                 Show this help message

Examples:
  ${SCRIPT_NAME} protein.pdb
  ${SCRIPT_NAME} protein.pdb -o mesh.vtk -v 0.5 --filter protein-only
  ${SCRIPT_NAME} --config config/protein_only.conf
  ${SCRIPT_NAME} --generate-config > config/default.conf
  ${SCRIPT_NAME} --batch file1.pdb,file2.pdb -p 1.5
EOF
}

die() { echo "Error: $1" >&2; exit 1; }

trim() {
    local s="$1"
    s="${s#"${s%%[![:space:]]*}"}"
    s="${s%"${s##*[![:space:]]}"}"
    echo "$s"
}

declare -A CONFIG_MAP
parse_config_file() {
    local file="$1"
    [[ -f "$file" ]] || die "Config file not found: $file"
    local section=""
    while IFS= read -r line || [[ -n "$line" ]]; do
        line="${line%%#*}"
        line="${line%%;*}"
        line=$(trim "$line")
        [[ -z "$line" ]] && continue
        if [[ "$line" =~ ^\[(.+)\]$ ]]; then
            section="${BASH_REMATCH[1]}"
            continue
        fi
        if [[ "$line" =~ ^([^=]+)=(.*)$ ]]; then
            local key
            key=$(trim "${BASH_REMATCH[1]}")
            local value
            value=$(trim "${BASH_REMATCH[2]}")
            CONFIG_MAP["${section}.${key}"]="$value"
        fi
    done <"$file"
}

set_from_config() {
    local key="$1" var_name="$2"
    if [[ -n "${CONFIG_MAP[$key]:-}" ]]; then
        printf -v "$var_name" "%s" "${CONFIG_MAP[$key]}"
    fi
}

log() { $quiet || echo "$@"; }
logv() { $verbose && ! $quiet && echo "$@"; }

generate_default_config() {
    cat <<'EOF'
[input]
file = example.pdb

[filter]
type = protein-only
keep_proteins = true
keep_water = false
keep_nucleic_acids = false
keep_ions = false
keep_lipids = true
keep_ligands = true

[voxelization]
voxel_size = 1.0
padding = 2.0

[output]
file = mesh.vtk
format = vtk

[options]
verbose = true
EOF
}

is_number() { [[ "$1" =~ ^-?[0-9]+([.][0-9]+)?$ ]]; }

allowed_filter() {
    case "$1" in
    none|all|protein-only|no-water|custom) return 0 ;;
    *) return 1 ;;
    esac
}

allowed_format() {
    case "$1" in
    vtk|gid|json|txt) return 0 ;;
    *) return 1 ;;
    esac
}

normalize_bool() {
    case "$(echo "$1" | tr '[:upper:]' '[:lower:]')" in
    1|true|yes|on) echo "true" ;;
    0|false|no|off|"") echo "false" ;;
    *) echo "$1" ;;
    esac
}

find_core_executable() {
    local build_dir="${SCRIPT_DIR}/build"
    local candidates=("biomesh2_core" "filter_workflow_example" "biomesh2_example" "filter_demo")
    for bin in "${candidates[@]}"; do
        if [[ -x "${build_dir}/${bin}" ]]; then
            CORE_CMD="${build_dir}/${bin}"
            CORE_NAME="$bin"
            return 0
        fi
    done
    die "No BioMesh2 executable found. Build the project in ${build_dir} first."
}

run_core_command() {
    local input_file="$1"
    local args=("$CORE_CMD")
    case "$CORE_NAME" in
    biomesh2_example)
        args+=("$input_file" "$padding")
        ;;
    filter_workflow_example|filter_demo)
        args+=("$input_file")
        ;;
    *)
        args+=("$input_file" "$padding")
        ;;
    esac
    logv "Running: ${args[*]}"
    if "${args[@]}"; then
        log "✓ Completed: $input_file"
    else
        die "Run failed for: $input_file"
    fi
}

# Defaults
input_file=""
output_file="mesh.vtk"
voxel_size="1.0"
padding="0.0"
filter_type="none"
output_format="vtk"
config_file=""
verbose=false
quiet=false
batch_inputs=()
input_from_config=false

# First pass: capture config file without mutating arguments
orig_args=("$@")
for ((i = 0; i < ${#orig_args[@]}; i++)); do
    case "${orig_args[$i]}" in
    -c|--config)
        if (( i + 1 < ${#orig_args[@]} )); then
            config_file="${orig_args[$((i + 1))]}"
        fi
        ;;
    --config=*)
        config_file="${orig_args[$i]#*=}"
        ;;
    esac
done

if [[ -z "$config_file" && -f "$DEFAULT_CONFIG" ]]; then
    config_file="$DEFAULT_CONFIG"
fi

if [[ -n "$config_file" ]]; then
    parse_config_file "$config_file"
    set_from_config "input.file" input_file
    set_from_config "filter.type" filter_type
    set_from_config "voxelization.voxel_size" voxel_size
    set_from_config "voxelization.padding" padding
    set_from_config "output.file" output_file
    set_from_config "output.format" output_format
    set_from_config "options.verbose" verbose
    set_from_config "options.quiet" quiet
    verbose=$(normalize_bool "$verbose")
    quiet=$(normalize_bool "$quiet")
    [[ -n "$input_file" ]] && input_from_config=true
fi

args=("$@")
idx=0
while [[ $idx -lt ${#args[@]} ]]; do
    arg="${args[$idx]}"
    case "$arg" in
    -i|--input)
        ((idx++)); input_file="${args[$idx]}"; input_from_config=false
        ;;
    --input=*)
        input_file="${arg#*=}"; input_from_config=false
        ;;
    -o|--output)
        ((idx++)); output_file="${args[$idx]}"
        ;;
    --output=*)
        output_file="${arg#*=}"
        ;;
    -v|--voxel-size)
        ((idx++)); voxel_size="${args[$idx]}"
        ;;
    --voxel-size=*)
        voxel_size="${arg#*=}"
        ;;
    -p|--padding)
        ((idx++)); padding="${args[$idx]}"
        ;;
    --padding=*)
        padding="${arg#*=}"
        ;;
    -f|--format)
        ((idx++)); output_format="${args[$idx]}"
        ;;
    --format=*)
        output_format="${arg#*=}"
        ;;
    --filter)
        ((idx++)); filter_type="${args[$idx]}"
        ;;
    --filter=*)
        filter_type="${arg#*=}"
        ;;
    -c|--config)
        ((idx++)) # already handled; skip value
        ;;
    --config=*)
        ;;
    --generate-config)
        generate_default_config
        exit 0
        ;;
    --batch)
        ((idx++)); batch_inputs+=("${args[$idx]}")
        ;;
    --batch=*)
        batch_inputs+=("${arg#*=}")
        ;;
    -V|--verbose)
        verbose=true
        ;;
    -q|--quiet)
        quiet=true
        ;;
    --version)
        echo "${SCRIPT_NAME} version ${VERSION}"
        exit 0
        ;;
    -h|--help)
        usage
        exit 0
        ;;
    --)
        ((idx++))
        for ((; idx < ${#args[@]}; idx++)); do
            batch_inputs+=("${args[$idx]}")
        done
        break
        ;;
    -*)
        die "Unknown option: $arg"
        ;;
    *)
        if $input_from_config || [[ -z "$input_file" ]]; then
            input_file="$arg"
            input_from_config=false
        else
            batch_inputs+=("$arg")
        fi
        ;;
    esac
    ((idx++))
done

expanded_batch=()
for item in "${batch_inputs[@]}"; do
    IFS=',' read -r -a parts <<<"$item"
    for part in "${parts[@]}"; do
        part=$(trim "$part")
        [[ -z "$part" ]] && continue
        expanded_batch+=("$part")
    done
done
batch_inputs=("${expanded_batch[@]}")

all_inputs=()
if [[ -n "$input_file" ]]; then
    all_inputs+=("$input_file")
fi
if [[ ${#batch_inputs[@]} -gt 0 ]]; then
    all_inputs+=("${batch_inputs[@]}")
fi

if [[ ${#all_inputs[@]} -eq 0 ]]; then
    usage
    die "No input PDB file specified."
fi

for f in "${all_inputs[@]}"; do
    [[ -f "$f" ]] || die "Input file not found: $f"
done

is_number "$voxel_size" || die "Voxel size must be numeric."
is_number "$padding" || die "Padding must be numeric."
allowed_filter "$filter_type" || die "Invalid filter type: $filter_type"
allowed_format "$output_format" || die "Invalid output format: $output_format"

print_banner

log "BioMesh2 configuration:"
log "  Inputs      : ${all_inputs[*]}"
log "  Output file : $output_file"
log "  Voxel size  : $voxel_size Å"
log "  Padding     : $padding Å"
log "  Filter type : $filter_type"
log "  Format      : $output_format"
log "  Verbose     : $verbose"
log "  Quiet       : $quiet"
log "  Config file : ${config_file:-<none>}"

find_core_executable
log "Using executable: $CORE_CMD"

for pdb in "${all_inputs[@]}"; do
    run_core_command "$pdb"
done

log "All tasks completed successfully."
