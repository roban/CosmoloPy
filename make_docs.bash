#!/bin/bash

epydoc --exclude='cosmolopy.EH._power' --exclude='cosmolopy.EH.power' -v --no-private --no-frames --html --docformat restructuredtext cosmolopy/ -o www/docAPI/