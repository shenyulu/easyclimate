Remove-Item .\test\baseline_images\*.png
pytest --mpl-generate-path="test/baseline_images"
