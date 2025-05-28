#!/bin/bash

# Duration in seconds (20 hours)
DURATION=$((20 * 60 * 60))

# Keep the computer awake
echo "Keeping the computer awake for 20 hours..."
caffeinate -u -t "$DURATION" &

# Wait for the specified duration
wait

echo "Done. The computer was kept awake for 20 hours."
