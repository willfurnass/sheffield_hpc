#! /bin/bash

# Outputs the number of GPUs requested by an SGE job.
function getSGERequestedGPUCount {
    export PATH="$SGE_ROOT/bin/lx-amd64/:$PATH"

    # Find job id even when using qrsh
    JOB_ID=${JOB_ID-$(cat /proc/$$/cpuset | cut -d/ -f3 | cut -d. -f1)}

    ngpus_per_core=$(SGE_SINGLE_LINE=1 qstat -j $JOB_ID | grep -Po '^hard resource_list.*gpu=\K[1-9][0-9]*')
    if [ -z $ngpus_per_core ]; then
        ngpus_per_core=0
    fi

    ncpus=1
    if $(qstat -j $JOB_ID | grep -Eq '^parallel environment: .*(smp|openmp).*'); then
        ncpus=$(qstat -j $JOB_ID | grep -Po 'parallel environment: .*(smp|openmp).* range: \K[1-9][0-9]*')
    fi

    ngpus=$(( $ncpus * $ngpus_per_core ))

    echo "$ngpus"
}

# Function to set the value of CUDA_VISIBLE_DEVICES to the id of the least-utilised device in the system.
# This will not set the device if CUDA_VISIBLE_DEVICES has already been set, to prevent issues on the ShARC DGX-1 (Node 126) or other private nodes which use the device allocation model
function setCUDAVisibileDevices {

    # If CUDA_VISIBLE_DEVICES is not currently set then proceed.
    if [ -z "$CUDA_VISIBLE_DEVICES" ]; then

        # If the nvidia-smi command is available
        if [ $( command -v nvidia-smi ) ]; then

            # Initialise the requested number of GPUs to 1
            NUM_GPUS_REQUESTED=1

            # If there is at least one argument, use the first argument
            if [ $# -gt 0 ]; then
                NUM_GPUS_REQUESTED=$1
            fi

            # If the value is not a positive integer, error.
            if [[ ! $NUM_GPUS_REQUESTED =~ ^[1-9][0-9]*$ ]]; then
                >&2 echo "Error: ${FUNCNAME[0]} only accepts integer values larger than 0"
                exit 1
            fi

            # Query nvidia-smi for device utilisation information
            nvidiasmi_output=$(nvidia-smi --query-gpu=index,name,memory.free,utilization.gpu,utilization.memory,temperature.gpu --format=csv,noheader,nounits)

            # Find the number of GPUs in this system/node
            devices_found=$(echo "$nvidiasmi_output" | wc -l )

            # If fewer devices are found than requested, error.
            if [ "$devices_found" -lt "$NUM_GPUS_REQUESTED" ]; then
                >&2 echo "Error: Requested $NUM_GPUS_REQUESTED but only found $devices_found"
                exit 1
            fi

            # Sort and filter the output of nvidia-smi to rank GPUs in reverse order of busyness.
            # Return a comma separated list of NUM_GPUS_REQUESTED device indices
            deviceid=$(echo "$nvidiasmi_output" | sort -t, -k3,3nr -k4,4n -k5,5n -k6,6n -k1,1nr | head -n $NUM_GPUS_REQUESTED | cut -d, -f1 | paste -sd "," -)

            # Output and export the cuda visible devices variable
            echo "Setting CUDA_VISIBLE_DEVICES=$deviceid"
            export CUDA_VISIBLE_DEVICES=$deviceid

            # Also export export CUDA_DEVICE_ORDER=PCI_BUS_ID, so that the cuda runtime enumerates in the same order as nvidia-smi, rather than in performance order. Requires CUDA 7 runtime or newer.
            export CUDA_DEVICE_ORDER=PCI_BUS_ID

        else
            >&2 echo "nvidia-smi is not available. Aborting."
        fi

    else
        # If CUDA_VISIBLE_DEVICES is already set, output the value to the user.
        >&2 echo "CUDA_VISIBLE_DEVICES is already set. CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES"
    fi
}


# No argument, defaults to 1 device
# --------
# setCUDAVisibileDevices

# Accepts 1 argument
# --------
# setCUDAVisibileDevices 2


# Can be used in conjunction with SGE to get the same number of devices as requested by the job.
# --------
setCUDAVisibileDevices $(getSGERequestedGPUCount)
