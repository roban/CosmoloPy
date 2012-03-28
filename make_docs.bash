#!/bin/bash

epydoc --exclude='cosmolopy.EH._power' --exclude='cosmolopy.EH._tf_fit' --exclude='cosmolopy.EH.power' --exclude='cosmolopy.EH.tf_fit' -v --no-private --no-frames --html --docformat restructuredtext cosmolopy/ -o www/docAPI/