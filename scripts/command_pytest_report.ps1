Remove-Item -Recurse -Force .\htmlcov
pytest --cov=easyclimate --cov-report=html
