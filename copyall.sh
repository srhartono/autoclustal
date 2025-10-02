ver=$(cat version.txt);
echo $ver
sed -i.bak -E "s/^version\s*=.*/version = \"$ver\"/" pyproject.toml

mkdir -p autoclustal
cp -fr ./bin ./autoclustal/
cp -fr ./modules ./autoclustal/
cp -fr ./conf ./autoclustal/
cp -f ./version.txt ./autoclustal/
cp -f ./pyproject.toml ./autoclustal/
cp -f ./README.md ./autoclustal/
cp -f ./environment.yml ./autoclustal/
cp -f ./nextflow.config ./autoclustal/

rm -fr ./dist/
rm -fr ./build/
rm -fr ./autoclustal.egg-info/
python -m build
python -m twine check dist/*
#python -m twine upload --repository testpypi dist/*
#python -m twine upload --repository pypi dist/*
#pip install -U autoclustal
echo "Build complete. You can now upload the package to PyPI using twine."