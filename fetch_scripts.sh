#!/bin/bash

REPO_URL_SSH="git@github.com:rcgsheffield/hpc-examples.git"
REPO_URL_HTTPS="https://github.com/rcgsheffield/hpc-examples.git"
REPO_DIR="hpc-examples"
USE_SSH=true  # Default to SSH unless --https is provided

# Check if --https was passed
if [[ "$1" == "--https" ]]; then
    USE_SSH=false
fi

# Function to check SSH access
check_ssh_access() {
    echo "Checking SSH connection to GitHub..."
    if ssh -o BatchMode=yes -T git@github.com 2>&1 | grep -q "successfully authenticated"; then
        echo "✅ SSH authentication successful."
        return 0
    else
        echo -e "\n⚠ SSH authentication to GitHub failed."
        echo "Possible reasons:"
        echo "- Your SSH key may not be added to GitHub."
        echo "- The SSH agent may not be running."
        echo "- Your firewall or network may block SSH connections."
        return 1
    fi
}

# Offer HTTPS fallback if SSH fails
clone_with_https() {
    read -p "Would you like to try cloning with HTTPS instead? (y/n): " response
    if [[ "$response" =~ ^[Yy]$ ]]; then
        echo "Cloning using HTTPS..."
        if git clone "$REPO_URL_HTTPS"; then
            echo -e "\n✅ Repository successfully cloned using HTTPS."
            echo -e "\n⚠ Note: If you want to **make changes and push them back**, SSH access is required."
            echo "To switch to SSH later, run:"
            echo "  cd $REPO_DIR"
            echo "  git remote set-url origin $REPO_URL_SSH"
        else
            echo "❌ HTTPS clone failed as well. Check your network and Git settings."
            exit 1
        fi
    else
        echo "Aborting operation. Please configure SSH and try again."
        exit 1
    fi
}

# Clone or update the repository
if [ ! -d "$REPO_DIR" ]; then
    if $USE_SSH; then
        echo "Attempting to clone the hpc-examples repository using SSH..."
        if check_ssh_access && git clone "$REPO_URL_SSH"; then
            echo "✅ Repository successfully cloned using SSH."
        else
            echo "❌ SSH clone failed."
            clone_with_https
        fi
    else
        echo "Cloning using HTTPS..."
        git clone "$REPO_URL_HTTPS" && echo "✅ Repository successfully cloned using HTTPS."
    fi
else
    echo "Updating the hpc-examples repository..."
    cd "$REPO_DIR" || exit
    git pull origin main
    cd ..
fi

echo "✅ Repository setup is complete."

