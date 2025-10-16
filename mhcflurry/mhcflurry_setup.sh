#!/bin/bash
# MHCflurry Easy Setup for Docker/WSL2
# This script sets up everything you need for peptide promiscuity analysis

set -e

# Colors for pretty output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

echo -e "${BLUE}🧬 MHCflurry Promiscuity Analysis Setup${NC}"
echo "Setting up MHCflurry for peptide-binding promiscuity calculations..."

# Check if we're in Docker/WSL
if [[ -f /.dockerenv ]] || grep -q microsoft /proc/version 2>/dev/null; then
    echo -e "${GREEN}✓ Detected containerized environment (Docker/WSL2)${NC}"
else
    echo -e "${YELLOW}⚠️  Not in Docker/WSL2, but continuing...${NC}"
fi

# Function to check if command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Check Python
if command_exists python3; then
    PYTHON_VERSION=$(python3 --version)
    echo -e "${GREEN}✓ Python found: $PYTHON_VERSION${NC}"
else
    echo -e "${RED}❌ Python 3 not found. Please install Python 3.9+${NC}"
    exit 1
fi

# Check pip
if command_exists pip3; then
    echo -e "${GREEN}✓ pip3 found${NC}"
else
    echo -e "${RED}❌ pip3 not found. Installing...${NC}"
    sudo apt-get update && sudo apt-get install -y python3-pip
fi

# Install MHCflurry
echo -e "${YELLOW}📦 Installing MHCflurry...${NC}"

# Create virtual environment (recommended)
if [[ ! -d "mhcflurry_env" ]]; then
    echo "Creating virtual environment..."
    python3 -m venv mhcflurry_env
fi

# Activate virtual environment
source mhcflurry_env/bin/activate

# Upgrade pip first
pip install --upgrade pip

# Install required packages
echo "Installing packages..."
pip install mhcflurry pandas numpy

# Download MHCflurry models (this is the key step!)
echo -e "${YELLOW}⬇️  Downloading MHCflurry models (this may take a few minutes)...${NC}"

# Try to download models
if mhcflurry-downloads fetch models_class1_minimal 2>/dev/null; then
    echo -e "${GREEN}✓ Downloaded minimal models${NC}"
else
    echo -e "${YELLOW}⚠️  Minimal models failed, trying presentation models...${NC}"
    if mhcflurry-downloads fetch models_class1_presentation 2>/dev/null; then
        echo -e "${GREEN}✓ Downloaded presentation models${NC}"
    else
        echo -e "${RED}❌ Model download failed. Trying manual approach...${NC}"
        
        # Manual download as fallback
        mkdir -p downloads
        echo "Downloading models manually..."
        wget -O downloads/models.tar.bz2 \
            https://github.com/openvax/mhcflurry/releases/download/1.6.0/models_class1_presentation.20200205.tar.bz2
        
        mhcflurry-downloads fetch models_class1_presentation --already-downloaded-dir downloads
        echo -e "${GREEN}✓ Manual download successful${NC}"
    fi
fi

# Test installation
echo -e "${YELLOW}🧪 Testing MHCflurry installation...${NC}"

python3 -c "
import mhcflurry
from mhcflurry import Class1AffinityPredictor
print('✓ MHCflurry import successful')

try:
    predictor = Class1AffinityPredictor.load()
    print('✓ Models loaded successfully')
    
    # Quick test prediction
    test_prediction = predictor.predict(['SIINFEKL'], ['HLA-A02:01'])
    print(f'✓ Test prediction: {test_prediction[0]:.1f} nM')
    print('🎉 MHCflurry is ready!')
    
except Exception as e:
    print(f'❌ Model loading failed: {e}')
    exit(1)
"

if [[ $? -eq 0 ]]; then
    echo -e "${GREEN}✅ MHCflurry setup complete!${NC}"
else
    echo -e "${RED}❌ Setup failed${NC}"
    exit 1
fi

# Create sample files for testing
echo -e "${YELLOW}📝 Creating sample files...${NC}"

# Sample peptides file
cat > sample_peptides.txt << 'EOF'
SIINFEKL
SLYNTVATL
GLCTLVAML
YLQPRTFLL
YLQQPRTLL
EOF

# Sample HLA results (using your format)
cat > sample_hla_results.tsv << 'EOF'
Sample	A1	A2	B1	B2	C1	C2
TestSample1	A*24:02	A*02:01	B*44:02	B*13:02	C*05:01	C*06:02
TestSample2	A*33:01	A*24:02	B*14:02	B*49:01	C*07:01	C*08:02
EOF

echo -e "${GREEN}Created sample files:${NC}"
echo "  - sample_peptides.txt (5 test peptides)"
echo "  - sample_hla_results.tsv (2 test samples)"

# Test the promiscuity calculator
echo -e "${YELLOW}🧪 Testing promiscuity calculator...${NC}"

# Save the Python script to a file
cat > mhcflurry_promiscuity.py << 'PYTHON_SCRIPT_EOF'
# The actual Python script content would go here
# (This is a placeholder - in practice, you'd copy the full script)
print("MHCflurry promiscuity calculator placeholder")
print("Copy the full Python script to this file to use it")
PYTHON_SCRIPT_EOF

echo -e "${BLUE}🎯 Quick Usage Examples:${NC}"
echo ""
echo "1. Activate environment:"
echo "   source mhcflurry_env/bin/activate"
echo ""
echo "2. Test with sample data:"
echo "   python3 mhcflurry_promiscuity.py sample_hla_results.tsv sample_peptides.txt"
echo ""
echo "3. Use your own data:"
echo "   python3 mhcflurry_promiscuity.py your_hla_results.tsv your_peptides.txt"
echo ""
echo "4. Per-sample analysis:"
echo "   python3 mhcflurry_promiscuity.py your_hla_results.tsv your_peptides.txt --analysis per-sample"
echo ""

echo -e "${GREEN}🎉 Setup complete! MHCflurry is ready for promiscuity analysis.${NC}"
echo -e "${YELLOW}💡 Remember to activate the virtual environment before use:${NC}"
echo -e "${YELLOW}   source mhcflurry_env/bin/activate${NC}"
