################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/MvBiClus.cpp \
../src/MvClus.cpp \
../src/MvLrmaL0.cpp \
../src/MvLrmaL1.cpp \
../src/MvSsvd.cpp \
../src/RcppExports.cpp \
../src/cluster.cpp \
../src/main.cpp \
../src/utils.cpp 

OBJS += \
./src/MvBiClus.o \
./src/MvClus.o \
./src/MvLrmaL0.o \
./src/MvLrmaL1.o \
./src/MvSsvd.o \
./src/RcppExports.o \
./src/cluster.o \
./src/main.o \
./src/utils.o 

CPP_DEPS += \
./src/MvBiClus.d \
./src/MvClus.d \
./src/MvLrmaL0.d \
./src/MvLrmaL1.d \
./src/MvSsvd.d \
./src/RcppExports.d \
./src/cluster.d \
./src/main.d \
./src/utils.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -std=c++0x -DINSIDE -I/usr/share/R/include -I/usr/lib/R/site-library/RcppArmadillo/include -I/usr/local/lib/R/site-library/Rcpp/include -I/usr/local/lib/R/site-library/Rcpp/include/Rcpp -I/usr/local/lib/R/site-library/RInside/include -O0 -g3 -Wall -c -fmessage-length=0 -m64 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


