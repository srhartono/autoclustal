$version = Get-Content .\version.txt | ForEach-Object { $_.Trim() }
(Get-Content pyproject.toml) -replace '^version\s*=.*', "version = `"$version`"" | Set-Content pyproject.toml
Write-Output "Updated pyproject.toml to version $version"

# 1. Make sure the folder exists (-Force makes parent dirs if needed)
New-Item -ItemType Directory -Force -Path .\autoclustal

# 2. Copy bin → autoclustal\bin
Copy-Item -Recurse -Force .\bin\* .\autoclustal\bin\

# 3. Copy modules → autoclustal\modules
Copy-Item -Recurse -Force .\modules\* .\autoclustal\modules\

# 4. Copy conf → autoclustal\conf
Copy-Item -Recurse -Force .\conf\* .\autoclustal\conf\

# 5. Copy README.md → autoclustal\
Copy-Item -Force .\pyproject.toml .\autoclustal\

# 5. Copy README.md → autoclustal\
Copy-Item -Force .\version.txt .\autoclustal\

# 5. Copy README.md → autoclustal\
Copy-Item -Force LICENSE .\autoclustal\
Copy-Item -Force MANIFEST.in .\autoclustal\

# 5. Copy README.md → autoclustal\
Copy-Item -Force .\README.md .\autoclustal\

# 6. Copy environment.yml → autoclustal\
Copy-Item -Force .\environment.yml .\autoclustal\

# 7. Copy nextflow.config → autoclustal\
Copy-Item -Force .\nextflow.config .\autoclustal\

if (Test-Path dist) { Remove-Item -Recurse -Force dist } ; if (Test-Path build) { Remove-Item -Recurse -Force build } ; if (Test-Path autoclustal.egg-info) { Remove-Item -Recurse -Force autoclustal.egg-info }
python -m build
python -m twine check dist/*
pip install .\dist\autoclustal-$version.tar.gz 
#python -m twine upload --repository testpypi dist/*
#python -m twine upload --repository pypi dist/*
# pip install -U autoclustal
Write-Output "Build complete. You can now upload the package to PyPI using twine."