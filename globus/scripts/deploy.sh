#!/bin/bash
# OneRoof Globus Action Provider Deployment Script

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
GLOBUS_DIR="$PROJECT_ROOT/globus"

echo "OneRoof Globus Action Provider Deployment"
echo "========================================="

# Check for required files
if [ ! -f "$GLOBUS_DIR/config/.env" ]; then
    echo "Error: config/.env not found!"
    echo "Please copy config/.env.template to config/.env and configure it."
    exit 1
fi

# Load environment variables
export $(grep -v '^#' "$GLOBUS_DIR/config/.env" | xargs)

# Validate required variables
required_vars=(
    "ONEROOF_PATH"
    "ONEROOF_WORK_DIR"
    "NEXTFLOW_PATH"
)

for var in "${required_vars[@]}"; do
    if [ -z "${!var:-}" ]; then
        echo "Error: $var is not set in config/.env"
        exit 1
    fi
done

# Create work directory if it doesn't exist
echo "Creating work directory: $ONEROOF_WORK_DIR"
mkdir -p "$ONEROOF_WORK_DIR"

# Install Python dependencies
echo "Installing Python dependencies..."
cd "$GLOBUS_DIR/action_provider"
if command -v pip3 &> /dev/null; then
    pip3 install -r requirements.txt
elif command -v pip &> /dev/null; then
    pip install -r requirements.txt
else
    echo "Error: pip not found. Please install Python 3 and pip."
    exit 1
fi

# Check if running as systemd service or standalone
if [ "${1:-}" == "--systemd" ]; then
    # Create systemd service
    echo "Creating systemd service..."

    SERVICE_FILE="/etc/systemd/system/oneroof-action-provider.service"

    if [ "$EUID" -ne 0 ]; then
        echo "Error: --systemd requires root privileges"
        echo "Run: sudo $0 --systemd"
        exit 1
    fi

    cat > "$SERVICE_FILE" << EOF
[Unit]
Description=OneRoof Globus Action Provider
After=network.target

[Service]
Type=simple
User=$SUDO_USER
WorkingDirectory=$GLOBUS_DIR/action_provider
Environment="PATH=/usr/local/bin:/usr/bin:/bin"
EnvironmentFile=$GLOBUS_DIR/config/.env
ExecStart=$(which python3) $GLOBUS_DIR/action_provider/oneroof_action_provider.py
Restart=always
RestartSec=10

[Install]
WantedBy=multi-user.target
EOF

    systemctl daemon-reload
    systemctl enable oneroof-action-provider
    systemctl start oneroof-action-provider

    echo "Service installed and started!"
    echo "Check status with: systemctl status oneroof-action-provider"
    echo "View logs with: journalctl -u oneroof-action-provider -f"

elif [ "${1:-}" == "--docker" ]; then
    # Deploy with Docker
    echo "Building Docker image..."

    cd "$GLOBUS_DIR"

    # Create Dockerfile if it doesn't exist
    if [ ! -f "Dockerfile" ]; then
        cat > Dockerfile << 'EOF'
FROM python:3.10-slim

WORKDIR /app

COPY action_provider/requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY action_provider/oneroof_action_provider.py .
COPY config/.env /app/.env

EXPOSE 5000

CMD ["python", "oneroof_action_provider.py"]
EOF
    fi

    docker build -t oneroof-action-provider .

    echo "Starting Docker container..."
    docker run -d \
        --name oneroof-action-provider \
        --restart unless-stopped \
        -p 5000:5000 \
        -v "$ONEROOF_PATH:$ONEROOF_PATH:ro" \
        -v "$ONEROOF_WORK_DIR:$ONEROOF_WORK_DIR" \
        --env-file "$GLOBUS_DIR/config/.env" \
        oneroof-action-provider

    echo "Container started!"
    echo "Check status with: docker ps"
    echo "View logs with: docker logs -f oneroof-action-provider"

else
    # Run in foreground
    echo "Starting action provider in foreground..."
    echo "Press Ctrl+C to stop"
    echo ""
    echo "For production deployment, use:"
    echo "  $0 --systemd    # Install as systemd service"
    echo "  $0 --docker     # Deploy with Docker"
    echo ""

    cd "$GLOBUS_DIR/action_provider"
    python3 oneroof_action_provider.py
fi
