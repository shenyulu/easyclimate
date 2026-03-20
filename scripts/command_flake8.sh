#!/bin/sh

# 获取当前 Python 解释器路径
# Gets the current Python interpreter path
PYTHON_PATH=$(which python)

# 检查是否找到 Python
# Check if Python is found
if [ -z "$PYTHON_PATH" ]; then
    echo "Error: Python not found in PATH."
    exit 1
fi

# 获取 Python 主次版本号（例如 "3.11"）
# Get Python primary and secondary version numbers (e.g. "3.11")
PYTHON_VERSION=$(python -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")

echo "Using Python: $PYTHON_PATH (Version: $PYTHON_VERSION)"
python -m flake8 src test scripts --max-line-length=88 --extend-ignore=E203,W503
