#!/bin/bash
# WSL2 Virtual Environment Setup for HLA Analysis Project

echo "Setting up Python virtual environment for HLA analysis..."

# Navigate to your project directory
cd /mnt/c/hla_analysis

# 1. Check if Python3 and pip are installed
echo "Checking Python installation..."
python3 --version
if [ $? -ne 0 ]; then
    echo "Python3 not found. Installing..."
    sudo apt update
    sudo apt install python3 python3-pip python3-venv -y
fi

# 2. Create virtual environment
echo "Creating virtual environment..."
python3 -m venv hla_analysis_env

# 3. Activate virtual environment
echo "Activating virtual environment..."
source hla_analysis_env/bin/activate

# 4. Upgrade pip
echo "Upgrading pip..."
pip install --upgrade pip

# 5. Install required packages
echo "Installing required Python packages..."
pip install pandas
pip install pathlib  # Usually included in Python 3.4+, but just in case
pip install openpyxl  # For Excel output if needed

# 6. Create a requirements.txt file for future reference
echo "Creating requirements.txt..."
cat > requirements.txt << EOF
pandas>=1.3.0
openpyxl>=3.0.0
pathlib-abc>=0.1.0
EOF

# 7. Install from requirements file (to verify)
pip install -r requirements.txt

echo ""
echo "Virtual environment setup complete!"
echo ""
echo "To activate the environment in the future, run:"
echo "  cd /mnt/c/hla_analysis"
echo "  source hla_analysis_env/bin/activate"
echo ""
echo "To deactivate when done:"
echo "  deactivate"
echo ""
echo "Your Python script should now run successfully!"
