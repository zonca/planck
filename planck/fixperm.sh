#!/bin/bash
python -m compileall .
chgrp -R cmb *
chmod -R g+wrX,o+rX *
