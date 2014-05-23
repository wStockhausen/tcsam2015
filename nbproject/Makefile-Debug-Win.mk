#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-MacOSX
CND_DLIB_EXT=dylib
CND_CONF=Debug-Win
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/src/CatchData.o \
	${OBJECTDIR}/src/ModelConfiguration.o \
	${OBJECTDIR}/src/ModelConstants.o \
	${OBJECTDIR}/src/ModelData.o \
	${OBJECTDIR}/src/ModelFunctions.o \
	${OBJECTDIR}/src/ModelIndexBlocks.o \
	${OBJECTDIR}/src/ModelIndexFunctions.o \
	${OBJECTDIR}/src/ModelParameterInfoTypes.o \
	${OBJECTDIR}/src/ModelParametersInfo.o \
	${OBJECTDIR}/src/ModelSelectivities.o \
	${OBJECTDIR}/src/SummaryFunctions.o \
	${OBJECTDIR}/src/TCSAM2015.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/tcsam2015

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/tcsam2015: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/tcsam2015 ${OBJECTFILES} ${LDLIBSOPTIONS}

${OBJECTDIR}/src/CatchData.o: src/CatchData.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/CatchData.o src/CatchData.cpp

${OBJECTDIR}/src/ModelConfiguration.o: src/ModelConfiguration.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/ModelConfiguration.o src/ModelConfiguration.cpp

${OBJECTDIR}/src/ModelConstants.o: src/ModelConstants.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/ModelConstants.o src/ModelConstants.cpp

${OBJECTDIR}/src/ModelData.o: src/ModelData.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/ModelData.o src/ModelData.cpp

${OBJECTDIR}/src/ModelFunctions.o: src/ModelFunctions.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/ModelFunctions.o src/ModelFunctions.cpp

${OBJECTDIR}/src/ModelIndexBlocks.o: src/ModelIndexBlocks.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/ModelIndexBlocks.o src/ModelIndexBlocks.cpp

${OBJECTDIR}/src/ModelIndexFunctions.o: src/ModelIndexFunctions.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/ModelIndexFunctions.o src/ModelIndexFunctions.cpp

${OBJECTDIR}/src/ModelParameterInfoTypes.o: src/ModelParameterInfoTypes.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/ModelParameterInfoTypes.o src/ModelParameterInfoTypes.cpp

${OBJECTDIR}/src/ModelParametersInfo.o: src/ModelParametersInfo.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/ModelParametersInfo.o src/ModelParametersInfo.cpp

${OBJECTDIR}/src/ModelSelectivities.o: src/ModelSelectivities.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/ModelSelectivities.o src/ModelSelectivities.cpp

${OBJECTDIR}/src/SummaryFunctions.o: src/SummaryFunctions.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/SummaryFunctions.o src/SummaryFunctions.cpp

${OBJECTDIR}/src/TCSAM2015.o: src/TCSAM2015.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/TCSAM2015.o src/TCSAM2015.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/tcsam2015

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
