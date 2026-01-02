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
  -p, --padding VALUE        Bounding box padding in Angstroms (default: 2.0)
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

trim() { echo "$1" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//'; }

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
# Placeholder input file; replace with your PDB path
file = input.pdb

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

append_batch_items() {
    local value="$1"
    IFS=',' read -r -a parts <<<"$value"
    for part in "${parts[@]}"; do
        part=$(trim "$part")
        [[ -z "$part" ]] && continue
        cli_batch+=("$part")
    done
}

find_core_executable() {
    local build_dir="${SCRIPT_DIR}/build"
    local candidates=("filter_workflow_example" "biomesh2_example" "filter_demo")
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
    local run_index="$2"
    local run_total="$3"
    local args=("$CORE_CMD")
    local env_vars=(
        "BIOMESH2_VOXEL_SIZE=$voxel_size"
        "BIOMESH2_PADDING=$padding"
        "BIOMESH2_FILTER=$filter_type"
        "BIOMESH2_OUTPUT_FORMAT=$output_format"
    )
    local effective_output="$output_file"
    local base="${output_file%.*}"
    local ext="$output_format"
    [[ "$output_file" == *.* ]] && ext="${output_file##*.}"
    if (( run_total > 1 )); then
        local stem
        stem=$(basename "${input_file%.*}")
        effective_output="${base}_${stem}_${run_index}.${ext}"
    fi
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
    if env "${env_vars[@]}" "BIOMESH2_OUTPUT_FILE=$effective_output" "${args[@]}"; then
        log "✓ Completed: $input_file"
        if [[ -n "$effective_output" ]]; then
            local summary_file="${effective_output}.summary"
            mkdir -p "$(dirname "$summary_file")" 2>/dev/null || true
            cat >"$summary_file" <<EOF
# BioMesh2 run summary
Input file: $input_file
Output file: $effective_output
Output format: $output_format
Voxel size: $voxel_size
Padding: $padding
Filter: $filter_type
Executable: $CORE_CMD
EOF
            logv "Wrote summary to $summary_file"
        fi
    else
        die "Run failed for: $input_file"
    fi
}

# Defaults
input_file=""
output_file="mesh.vtk"
voxel_size="1.0"
padding="2.0"
filter_type="none"
output_format="vtk"
config_file=""
verbose=false
quiet=false
batch_inputs=()
cli_input=""
cli_output=""
cli_voxel=""
cli_padding=""
cli_filter=""
cli_format=""
cli_batch=()
cli_verbose_set=false
cli_quiet_set=false

[[ -f "$DEFAULT_CONFIG" ]] && config_file="$DEFAULT_CONFIG"

require_value() {
    local option="$1"
    local next_index="$2"
    local total="$3"
    (( next_index < total )) || die "Option ${option} requires a value."
}

args=("$@")
idx=0
while [[ $idx -lt ${#args[@]} ]]; do
    arg="${args[$idx]}"
    case "$arg" in
    -i|--input)
        require_value "$arg" $((idx + 1)) ${#args[@]}
        ((idx++)); cli_input="${args[$idx]}"
        ;;
    --input=*)
        cli_input="${arg#*=}"
        ;;
    -o|--output)
        require_value "$arg" $((idx + 1)) ${#args[@]}
        ((idx++)); cli_output="${args[$idx]}"
        ;;
    --output=*)
        cli_output="${arg#*=}"
        ;;
    -v|--voxel-size)
        require_value "$arg" $((idx + 1)) ${#args[@]}
        ((idx++)); cli_voxel="${args[$idx]}"
        ;;
    --voxel-size=*)
        cli_voxel="${arg#*=}"
        ;;
    -p|--padding)
        require_value "$arg" $((idx + 1)) ${#args[@]}
        ((idx++)); cli_padding="${args[$idx]}"
        ;;
    --padding=*)
        cli_padding="${arg#*=}"
        ;;
    -f|--format)
        require_value "$arg" $((idx + 1)) ${#args[@]}
        ((idx++)); cli_format="${args[$idx]}"
        ;;
    --format=*)
        cli_format="${arg#*=}"
        ;;
    --filter)
        require_value "$arg" $((idx + 1)) ${#args[@]}
        ((idx++)); cli_filter="${args[$idx]}"
        ;;
    --filter=*)
        cli_filter="${arg#*=}"
        ;;
    -c|--config)
        require_value "$arg" $((idx + 1)) ${#args[@]}
        ((idx++)); config_file="${args[$idx]}"
        ;;
    --config=*)
        config_file="${arg#*=}"
        ;;
    --generate-config)
        generate_default_config
        exit 0
        ;;
    --batch)
        require_value "$arg" $((idx + 1)) ${#args[@]}
        ((idx++)); append_batch_items "${args[$idx]}"
        ;;
    --batch=*)
        append_batch_items "${arg#*=}"
        ;;
    -V|--verbose)
        cli_verbose_set=true
        verbose=true
        ;;
    -q|--quiet)
        cli_quiet_set=true
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
            cli_batch+=("${args[$idx]}")
        done
        break
        ;;
    -*)
        die "Unknown option: $arg"
        ;;
    *)
        if [[ -z "$cli_input" ]]; then
            cli_input="$arg"
        else
            append_batch_items "$arg"
        fi
        ;;
    esac
    ((idx++))
done

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
fi

[[ -n "$cli_input" ]] && input_file="$cli_input"
[[ ${#cli_batch[@]} -gt 0 ]] && batch_inputs=("${cli_batch[@]}")
[[ -n "$cli_output" ]] && output_file="$cli_output"
[[ -n "$cli_voxel" ]] && voxel_size="$cli_voxel"
[[ -n "$cli_padding" ]] && padding="$cli_padding"
[[ -n "$cli_filter" ]] && filter_type="$cli_filter"
[[ -n "$cli_format" ]] && output_format="$cli_format"
$cli_verbose_set && verbose=true
$cli_quiet_set && quiet=true

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
log "  Voxel size  : $voxel_size Angstroms"
log "  Padding     : $padding Angstroms"
log "  Filter type : $filter_type"
log "  Format      : $output_format"
log "  Verbose     : $verbose"
log "  Quiet       : $quiet"
log "  Config file : ${config_file:-<none>}"

find_core_executable
log "Using executable: $CORE_CMD"

for idx in "${!all_inputs[@]}"; do
    run_core_command "${all_inputs[$idx]}" "$((idx + 1))" "${#all_inputs[@]}"
done

log "All tasks completed successfully."
