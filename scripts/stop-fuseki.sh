#!/bin/bash
# Stop Apache Jena Fuseki server

PID_FILE="/tmp/fuseki-gene.pid"

if [ -f "$PID_FILE" ]; then
    PID=$(cat "$PID_FILE")
    if kill -0 "$PID" 2>/dev/null; then
        echo "Stopping Fuseki (PID: $PID)..."
        kill "$PID"
        rm -f "$PID_FILE"
        echo "Stopped."
    else
        echo "Fuseki not running (stale PID file)"
        rm -f "$PID_FILE"
    fi
else
    # Try to find and kill any running fuseki
    PIDS=$(pgrep -f "fuseki-server.jar")
    if [ -n "$PIDS" ]; then
        echo "Found Fuseki processes: $PIDS"
        kill $PIDS
        echo "Stopped."
    else
        echo "Fuseki not running"
    fi
fi
