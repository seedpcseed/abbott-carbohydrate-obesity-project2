#!/bin/bash
# watch-build.sh - Abbott-Lurie ZC07 Next Steps Proposal

# Always run from this script's directory so relative paths work
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR" || exit 1

SOURCE_FILE="docs/zc07-next-steps-proposal.2026-02-01.md"
SOURCE_DIR="docs"
SOURCE_BASENAME="zc07-next-steps-proposal.2026-02-01.md"
OUTPUT="render/zc07-next-steps.pdf"

# Check if inotifywait is available
if ! command -v inotifywait &> /dev/null; then
    echo "Error: inotifywait is not installed."
    echo "Please install it with: sudo apt-get install inotify-tools"
    exit 1
fi

echo "Watching for changes to:"
echo "  - $SOURCE_FILE"
echo "Output: $OUTPUT"
echo ""

# Verify source file exists before watching
if [[ ! -f "$SOURCE_FILE" ]]; then
    echo "Error: Source file not found: $SOURCE_FILE"
    exit 1
fi
echo "Source file verified: $SOURCE_FILE"

# Debounce delay in seconds (adjust as needed)
DEBOUNCE_DELAY=${DEBOUNCE_DELAY:-3}
PID_FILE="/tmp/watch-build-zc07-$(basename $(pwd)).pid"

echo "Debounce delay: ${DEBOUNCE_DELAY}s (set DEBOUNCE_DELAY environment variable to change)"

# Check for PDF refresh tools
REFRESH_TOOLS_AVAILABLE=0
if command -v xdotool >/dev/null 2>&1 || command -v gdbus >/dev/null 2>&1 || command -v qdbus >/dev/null 2>&1; then
    REFRESH_TOOLS_AVAILABLE=1
fi

if [[ "${NO_AUTO_REFRESH:-}" != "1" ]] && [[ "${NO_AUTO_REFRESH:-}" != "true" ]]; then
    if [[ $REFRESH_TOOLS_AVAILABLE -eq 1 ]] || grep -qi microsoft /proc/version 2>/dev/null; then
        echo "PDF auto-refresh: enabled"
        if [[ $REFRESH_TOOLS_AVAILABLE -eq 0 ]] && grep -qi microsoft /proc/version 2>/dev/null; then
            echo "  (Windows PDF viewer typically auto-refreshes on file change)"
        fi
    else
        echo "PDF auto-refresh: disabled (no refresh tools found)"
        echo "  Install xdotool for Zathura refresh: sudo apt-get install xdotool"
    fi
else
    echo "PDF auto-refresh: disabled (NO_AUTO_REFRESH set)"
fi
echo ""

# Cleanup function
cleanup() {
    if [[ -f "$PID_FILE" ]]; then
        PENDING_PID=$(cat "$PID_FILE" 2>/dev/null)
        if [[ -n "$PENDING_PID" ]] && kill -0 "$PENDING_PID" 2>/dev/null; then
            kill "$PENDING_PID" 2>/dev/null
            wait "$PENDING_PID" 2>/dev/null
        fi
        rm -f "$PID_FILE"
    fi
    # Close file descriptor if open
    exec 3<&- 2>/dev/null
}

# Trap cleanup on exit
trap cleanup EXIT INT TERM

# Function to refresh PDF viewer if open
refresh_pdf() {
    local pdf_path="$1"
    local abs_path=$(readlink -f "$pdf_path" 2>/dev/null || realpath "$pdf_path" 2>/dev/null || echo "$pdf_path")
    
    # Skip refresh if disabled via environment variable
    if [[ "${NO_AUTO_REFRESH:-}" == "1" ]] || [[ "${NO_AUTO_REFRESH:-}" == "true" ]]; then
        return 0
    fi
    
    # Only try GUI refresh if DISPLAY is set (X11/GUI environment)
    if [[ -z "$DISPLAY" ]] && ! grep -qi microsoft /proc/version 2>/dev/null; then
        return 0  # No GUI environment, skip
    fi
    
    # Try to refresh different PDF viewers
    # Zathura (supports reload via keybinding 'r' key)
    if pgrep -x zathura >/dev/null 2>&1 && command -v xdotool >/dev/null 2>&1; then
        # Find zathura window with this PDF and send reload command
        local zathura_window=$(xdotool search --class zathura 2>/dev/null | head -1)
        if [[ -n "$zathura_window" ]]; then
            xdotool windowactivate --sync "$zathura_window" key r 2>/dev/null && echo "  ↻ Refreshed Zathura viewer"
            return 0
        fi
    fi
    
    # Okular (can reload via dbus)
    if pgrep -x okular >/dev/null 2>&1 && command -v qdbus >/dev/null 2>&1; then
        # Try to reload via dbus
        qdbus org.kde.okular-* /okularMainWindow_1 reload 2>/dev/null && echo "  ↻ Refreshed Okular viewer"
        return 0
    fi
    
    # Evince/GNOME Document Viewer (can reload via dbus)
    if pgrep -x evince >/dev/null 2>&1 && command -v gdbus >/dev/null 2>&1; then
        # Try to reload evince
        gdbus call --session --dest org.gnome.Evince --object-path /org/gnome/Evince/Window/0 --method org.gtk.Actions.Activate 'reload' [] {} 2>/dev/null && echo "  ↻ Refreshed Evince viewer"
        return 0
    fi
    
    # MuPDF (glfw viewer - can't easily refresh, but we can note it)
    if pgrep -f "mupdf.*$(basename "$pdf_path")" >/dev/null 2>&1; then
        echo "  ℹ MuPDF viewer detected - press 'r' to reload manually"
        return 0
    fi
    
    # Windows PDF viewer (if using WSL with Windows explorer)
    if grep -qi microsoft /proc/version 2>/dev/null; then
        # In WSL, PDFs opened via Windows viewer might auto-refresh, but we can touch the file
        touch "$pdf_path" 2>/dev/null
        # Note: Windows PDF viewers typically auto-refresh when file changes
        # Some PDF viewers in WSL might need manual refresh
        echo "  ℹ PDF updated - Windows viewer may auto-refresh"
    fi
}

# Build function
build_document() {
    local changed_file="$1"
    echo "[$(date +%H:%M:%S)] Building (triggered by $changed_file)..."
    
    # Capture pandoc output and exit code
    # Check if highlight.lua exists, use it if available
    LUA_FILTER_ARGS=()
    if [[ -f "${SCRIPT_DIR}/highlight.lua" ]]; then
        LUA_FILTER_ARGS=( "--lua-filter=${SCRIPT_DIR}/highlight.lua" )
    fi

    PANDOC_OUTPUT=$(pandoc "$SOURCE_FILE" -o "$OUTPUT" --pdf-engine=xelatex "${LUA_FILTER_ARGS[@]}" 2>&1)
    PANDOC_EXIT=$?
    
    # Show errors/warnings if any
    if echo "$PANDOC_OUTPUT" | grep -qi "error\|warning"; then
        echo "$PANDOC_OUTPUT" | grep -i "error\|warning"
    fi
    
    # Check if pandoc succeeded AND output file exists
    if [[ $PANDOC_EXIT -eq 0 ]] && [[ -f "$OUTPUT" ]]; then
        echo "[$(date +%H:%M:%S)] ✓ Build successful - Output: $OUTPUT"
        # Attempt to refresh PDF viewer
        refresh_pdf "$OUTPUT"
    else
        echo "[$(date +%H:%M:%S)] ✗ Build failed (exit code: $PANDOC_EXIT)"
        if [[ ! -f "$OUTPUT" ]]; then
            echo "  Output file was not created: $OUTPUT"
        fi
        if [[ -n "$PANDOC_OUTPUT" ]]; then
            echo "  Pandoc output:"
            echo "$PANDOC_OUTPUT" | sed 's/^/    /'
        fi
    fi
    echo ""
    
    # Clean up PID file after build completes
    rm -f "$PID_FILE"
}

# Function to handle file changes with debouncing
handle_file_change() {
    local filename="$1"
    
    # Check if the changed file matches our source file
    if [[ "$filename" == "$SOURCE_BASENAME" ]]; then
        echo "[$(date +%H:%M:%S)] Change detected in $SOURCE_FILE, waiting ${DEBOUNCE_DELAY}s..."
        
        # Kill any pending build
        if [[ -f "$PID_FILE" ]]; then
            PENDING_PID=$(cat "$PID_FILE" 2>/dev/null)
            if [[ -n "$PENDING_PID" ]] && kill -0 "$PENDING_PID" 2>/dev/null; then
                kill "$PENDING_PID" 2>/dev/null
                wait "$PENDING_PID" 2>/dev/null
            fi
            rm -f "$PID_FILE"
        fi
        
        # Start a new delayed build in the background
        (
            sleep "$DEBOUNCE_DELAY"
            build_document "$SOURCE_FILE"
        ) &
        echo $! > "$PID_FILE"
    fi
}

# Create render directory if it doesn't exist
mkdir -p render

# Do initial build
echo "Performing initial build..."
build_document "$SOURCE_FILE"

echo "Watching $SOURCE_DIR/ for changes..."
echo ""

# Watch the docs directory for our source file
# Use exec to avoid subshell issues with pipes
exec 3< <(inotifywait -m -e close_write,moved_to --format '%f' "$SOURCE_DIR" 2>/dev/null)
while IFS= read -r filename <&3; do
    handle_file_change "$filename"
done
exec 3<&-
