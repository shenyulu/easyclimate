pytest --mpl --mpl-baseline-path="test/baseline_images" --cov src

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

# 设置环境变量，强制 pre-commit 使用当前 Python
# Set environment variables to force pre-commit to use the current Python
$env:PRE_COMMIT_USE_SYSTEM_PYTHON = "1"
$env:VIRTUALENV_PYTHON = "$PYTHON_PATH"

# 运行 pre-commit
# Run pre-commit
Write-Host "Using Python: $PYTHON_PATH (Version: $PYTHON_VERSION)"
pre-commit run --all-files

python -u ".\scripts\check_deprecations.py"
