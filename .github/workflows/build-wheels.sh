#!/bin/bash
set -e -x

# Install swig3
yum localinstall -y http://springdale.math.ias.edu/data/puias/computational/6/x86_64//swig3012-3.0.12-3.sdl6.x86_64.rpm
ln /usr/local/swig/3.0.12/bin/swig /usr/bin/swig --symbolic
ln /usr/local/swig/3.0.12/bin/ccache-swig /usr/bin/ccache-swig --symbolic

# Compile wheels
for PYBIN in /opt/python/*[23][678]*/bin; do
    "${PYBIN}/pip" install -r /io/requirements.txt
    "${PYBIN}/pip" wheel /io/ -w wheelhouse/
done

# Bundle external shared libraries into the wheels
for whl in wheelhouse/cosmolopy*.whl; do
    auditwheel repair "$whl" --plat $PLAT -w /io/wheelhouse/
done

# Install packages and test
for PYBIN in /opt/python/*[23][678]*/bin/; do
    "${PYBIN}/pip" install cosmolopy --no-index -f /io/wheelhouse
    "${PYBIN}/python" -c "import cosmolopy; import cosmolopy.EH.power"
done
