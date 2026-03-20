# 获取当前 Python 解释器路径
# Gets the current Python interpreter path
$PYTHON_PATH = (Get-Command python).Source

# 检查是否找到 Python
# Check if Python is found
if (-not $PYTHON_PATH) {
    Write-Host "Error: Python not found in PATH."
    exit 1
}

# 获取 Python 主次版本号（例如 "3.11"）
# Get Python primary and secondary version numbers (e.g. "3.11")
$PYTHON_VERSION = & python -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')"

Write-Host "Using Python: $PYTHON_PATH (Version: $PYTHON_VERSION)"
Write-Host "Checking duplicate imports with flake8 (F811)..."
python -m flake8 src test scripts --select=F811
