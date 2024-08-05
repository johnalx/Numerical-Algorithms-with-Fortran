@echo off

set gccload=30

if not "%_4VER%" == "" echo.|time >fnumzeit.aus
if     "%_4VER%" == "" time < nul >fnumzeit.aus

dmake >dmake.aus 2>&1

if not "%_4VER%" == "" echo.|time >>fnumzeit.aus
if     "%_4VER%" == "" time < nul >>fnumzeit.aus

dmake test <makefile.in >makefile.aus 2>&1

if not "%_4VER%" == "" echo.|time >>fnumzeit.aus
if     "%_4VER%" == "" time < nul >>fnumzeit.aus
