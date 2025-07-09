# OneRoof Globus Flows Integration

This directory contains the Globus Flows orchestration layer for the OneRoof bioinformatics pipeline, enabling automated, scalable execution of sequencing analysis workflows across Globus endpoints.

## Overview

The OneRoof Globus Flows adapter provides:

- **Automated data transfers** between Globus endpoints (input data, reference files, results)
- **Pipeline orchestration** for both Nanopore and Illumina sequencing platforms
- **Error handling and notifications** for failed runs
- **Workspace management** with automatic cleanup
- **Flexible configuration** through environment variables and dotfiles

## Architecture

```
globus/
├── flows/                    # Flow definitions
│   └── oneroof_flow.json    # Main workflow definition
├── action_provider/          # Custom action provider
│   ├── oneroof_action_provider.py
│   └── requirements.txt
├── config/                   # Configuration templates
│   ├── .env.template
│   └── action_provider_config.json
├── scripts/                  # Deployment utilities
│   ├── deploy.sh
│   ├── register_flow.py
│   └── test_flow.py
└── README.md
```

## Prerequisites

1. **Globus Account and Endpoints**:
   - Active Globus account with appropriate permissions
   - Source endpoint (where your data is stored)
   - Compute endpoint (where OneRoof will run)
   - Destination endpoint (where results will be stored)

2. **OneRoof Pipeline Requirements**:
   - Nextflow installed on compute endpoint
   - Docker/Singularity for containerized execution (or Pixi for containerless)
   - CUDA-capable GPU for Nanopore basecalling (if using pod5 input)

3. **Globus Flows Requirements**:
   - Globus Flows subscription
   - Globus CLI installed locally
   - Python 3.8+ for action provider

## Quick Start

### 1. Clone and Configure

```bash
# Clone the OneRoof repository to your compute endpoint
git clone https://github.com/nrminor/oneroof.git
cd oneroof/globus

# Copy environment template
cp config/.env.template config/.env

# Edit config/.env with your values
# IMPORTANT: Never commit .env to version control. If pre-commits are installed, there is
# a pre-commit hook that will reject any commits with files ending with `.env`
```

### 2. Configure Environment Variables

Edit `config/.env` with your specific values:

```bash
# Required Globus settings
GLOBUS_CLIENT_ID=your-client-id
GLOBUS_CLIENT_SECRET=your-client-secret
ONEROOF_ACTION_URL=https://your-compute-endpoint.example.com:5000
ONEROOF_ACTION_SCOPE=https://auth.globus.org/scopes/your-scope

# Endpoint IDs
COMPUTE_ENDPOINT_ID=your-compute-endpoint-uuid
STORAGE_ENDPOINT_ID=your-storage-endpoint-uuid

# OneRoof paths (**absolute paths required**)
ONEROOF_PATH=/absolute/path/to/oneroof
ONEROOF_WORK_DIR=/absolute/path/to/work/directory
```

### 3. Deploy Action Provider

On your compute endpoint:

```bash
cd globus/scripts
./deploy.sh

# Or manually with good ol' pip:
cd ../action_provider
pip install -r requirements.txt
python oneroof_action_provider.py

# uv is also supported:
cd ../action_provider
uv run oneroof_action_provider.py

```

### 4. Register the Flow

```bash
cd globus/scripts
python register_flow.py

# This will output your flow ID
# Save this ID for running flows
```

### 5. Run a Flow

```bash
# Using the test script
python test_flow.py --flow-id YOUR_FLOW_ID \
  --platform nanopore \
  --input-data /path/on/source/endpoint/pod5_files \
  --primer-bed /path/to/primers.bed \
  --refseq /path/to/reference.fasta \
  --ref-gbk /path/to/reference.gbk \
  --output-path /path/on/destination/endpoint/results

# Or using Globus CLI directly
globus flows run YOUR_FLOW_ID --input flow_input.json
```

## Flow Input Schema

The flow expects the following input parameters:

```jsonc
{
  "run_id": "unique_run_identifier",
  "platform": "nanopore",  // or "illumina"
  "input_data": "/path/on/source/endpoint/to/data",
  "primer_bed": "/path/to/primers.bed",
  "refseq": "/path/to/reference.fasta",
  "ref_gbk": "/path/to/reference.gbk",
  "source_endpoint_id": "uuid-of-source-endpoint",
  "compute_endpoint_id": "uuid-of-compute-endpoint",
  "destination_endpoint_id": "uuid-of-destination-endpoint",
  "output_path": "/path/on/destination/endpoint/for/results",
  "nextflow_params": {
    "min_variant_frequency": 0.1,
    "downsample_to": 200
  },
  "profile": "docker",  // or "singularity", "containerless"
  "notification_recipients": ["email@example.com"],
  "keep_logs": false
}
```

## Advanced Configuration

### Using Pixi Environment

If you prefer to use Pixi instead of containers:

```bash
# In config/.env
USE_PIXI=true
PIXI_PATH=/path/to/pixi

# Ensure pixi environment is set up on compute endpoint
cd /path/to/oneroof
pixi install --frozen
```

### GPU Configuration

For Nanopore basecalling with GPU:

```bash
# In config/.env
ENABLE_GPU=true

# Ensure CUDA is available on compute endpoint
# The action provider will pass appropriate flags to Nextflow
```

### Custom Nextflow Parameters

Pass additional parameters through the flow input:

```jsonc
{
  "nextflow_params": {
    "pod5_batch_size": 100,
    "model": "dna_r10.4.1_e8.2_400bps_sup@v4.3.0",
    "low_memory": true,
    "skip_sylph": true
  }
}
```

### Security Considerations

1. **Environment Files**: Never commit `.env` files to version control
2. **API Keys**: Use the optional API_KEY setting for additional security
3. **Endpoint Access**: Ensure proper Globus ACLs on all endpoints
4. **Network Security**: Use HTTPS for action provider endpoints

### Monitoring and Debugging

1. **Check Flow Status**:
   ```bash
   globus flows run-status FLOW_ID RUN_ID
   ```

2. **View Action Provider Logs**:
   ```bash
   # On compute endpoint
   tail -f /path/to/oneroof_work_dir/RUN_ID/logs/nextflow.log
   ```

3. **Debug Failed Runs**:
   - Check notification emails for error summaries
   - Review logs in work directory (if keep_logs=true)
   - Use the Globus web interface to inspect flow execution

## Troubleshooting

### Common Issues

1. **Action Provider Not Reachable**:
   - Verify firewall rules allow incoming connections
   - Check ONEROOF_ACTION_URL is publicly accessible
   - Ensure SSL certificates are valid (if using HTTPS)

2. **Pipeline Fails to Start**:
   - Verify ONEROOF_PATH points to repository root
   - Check Nextflow is in PATH on compute endpoint
   - Ensure work directory has sufficient permissions

3. **Transfer Failures**:
   - Verify endpoint IDs are correct
   - Check Globus endpoint activation status
   - Ensure paths exist on source/destination endpoints

4. **Memory/Resource Issues**:
   - Use `low_memory` flag in nextflow_params
   - Adjust `pod5_batch_size` for GPU memory
   - Set `downsample_to` to limit coverage depth

### Getting Help

1. Check OneRoof documentation: [main README](../README.md)
2. Review Globus Flows documentation: https://docs.globus.org/flows/
3. Check action provider logs for detailed error messages
4. Open an issue in the OneRoof repository with:
   - Flow run ID
   - Error messages from notifications
   - Relevant log excerpts

## Development

### Testing Locally

```bash
# Run action provider in debug mode
cd globus/action_provider
export FLASK_DEBUG=true
python oneroof_action_provider.py

# Test individual actions
curl -X POST http://localhost:5000 \
  -H "Content-Type: application/json" \
  -d '{"action": "validate", "platform": "nanopore", ...}'
```

### Extending the Flow

To add new capabilities:

1. Add new states to `flows/oneroof_flow.json`
2. Implement corresponding handlers in `action_provider/oneroof_action_provider.py`
3. Update input schema as needed
4. Test thoroughly before deploying

## License

This Globus integration follows the same license as the OneRoof project.
