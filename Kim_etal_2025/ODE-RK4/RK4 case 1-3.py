#!/usr/bin/env python3

import numpy as np
from pandas import DataFrame

# Input numbers
dm1=0.00303;		# co-transcriptional degradation rate (unit: 1/s)					# kd1 in the paper
dm2=0.00708;		# post-transcriptional degradation rate (unit: 1/s)					# kd2 in the paper
tt5=28.2;			# 5' transcription time (unit: s)									# T5' in the paper
tt3=210;			# 3' transcription time (unit: s)									# T3' in the paper
T1_5=128.8;			# when the last RNAP passes 5' probing region (unit: s)				# t5' in the paper
T2=228.1;			# when the first RNAP passes transcription terminator (unit: s)		# T2 in the supplementary
T3=328.7;	# when the last RNAP passes transcription terminator (unit: s)		# T3 in the supplementary
T1_3=T3-(T2-tt3);		# when the last RNAP passes 3' probing region (unit: s)				# t3' in the paper
km1=0.0176;		# mRNA synthesis rate with an inducer (unit: 1/s)                      #ki in the supplementary
km2=0.000107;		# synthesis rate of basal level (unit: 1/s)                        #kbasal in the supplementary

# (5') Nascent mRNA 1 (steady-state mRNA level of basal level expression)
def mn1_5(t):
	return (km2/dm1)*(1-(np.exp(-dm1*(T2-tt5))))

# (5') Released mRNA 1 (steady-state mRNA level of basal level expression)
def mr1_5(t):
	return (km2/dm2)*(np.exp(-dm1*(T2-tt5)))

# (3') Nascent mRNA 1 (steady-state mRNA level of basal level expression)
def mn1_3(t):
	return (km2/dm1)*(1-(np.exp(-dm1*(T2-tt3))))

# (3') Released mRNA 1 (steady-state mRNA level of basal level expression)
def mr1_3(t):
	return (km2/dm2)*(np.exp(-dm1*(T2-tt3)))

# ODE n1 (5') (tt5 < t < T1_5)
def un1_5(t, y):
	return km1 - y * dm1

# ODE n2 (5') (T1_5 < t < T2)
def un2_5(t, y):
	return 0 - y * dm1

# ODE n3 (5') (T2 < t < T3)
def un3_5(t, y):
	return -km1*np.exp(-dm1*(T2-tt5)) - y * dm1

# ODE n4 (5') (T3 < t)
def un4_5(t, y):
	return 0

# ODE r1 (5') (tt5 < t < T1_5)
def ur1_5(t, y):
	return 0

# ODE r2 (5') (T1_5 < t < T2)
def ur2_5(t, y):
	return 0

# ODE r3 (5') (T2 < t < T3)
def ur3_5(t, y):
	return km1*np.exp(-dm1*(T2-tt5)) - y * dm2

# ODE r4 (5') (T3 < t)
def ur4_5(t, y):
	return 0 - y * dm2

# ODE n1 (3') (tt3 < t < T2)
def un1_3(t, y):
	return km1 - y * dm1

# ODE n2 (3') (T2 < t < T1_3)
def un2_3(t, y):
	return 0

# ODE n3 (3') (T1_3 < t < T3)
def un3_3(t, y):
	return -km1*np.exp(-dm1*(T2-tt3)) - y * dm1

# ODE n4 (3') (T3 < t)
def un4_3(t, y):
	return 0

# ODE r1 (3') (tt3 < t < T2)
def ur1_3(t, y):
	return 0

# ODE r2 (3') (T2 < t < T1_3)
def ur2_3(t, y):
	return km1*np.exp(-dm1*(T2-tt3)) - y * dm2

# ODE r3 (3') (T1_3 < t < T3)
def ur3_3(t, y):
	return km1*np.exp(-dm1*(T2-tt3)) - y * dm2

# ODE r4 (3') (T3 < t)
def ur4_3(t, y):
	return - y * dm2

# Initial condition (5')
yn0_5 = mn1_5(tt5)
yr0_5 = mr1_5(tt5)

# Initial condition (3')
yn0_3 = mn1_3(tt3)
yr0_3 = mr1_3(tt3)

# RK4 (5' nascent mRNA)
def rk4n_5(tt5, T1_5, T2, T3, T_last, h):
	num_total = np.arange(0, T_last, h)
	n_total = len(num_total)
	num_0 = np.arange(0, tt5, h)
	n_0 = len(num_0)
	num_1 = np.arange(tt5, T1_5, h)
	n_1 = len(num_1)
	num_2 = np.arange(T1_5, T2, h)
	n_2 = len(num_2)
	num_3 = np.arange(T2, T3, h)
	n_3 = len(num_3)
	num_4 = np.arange(T3, T_last, h)
	n_4 = len(num_4)
	time1 = [0] * (n_total + 1)
	RK4n_5 = [0] * (n_total + 1)
	time1[0] = t = 0
	RK4n_5[0] = y = yn0_5
	for i0 in range(1, n_0 + 1):
		time1[i0] = t = time1[i0-1] + h
		RK4n_5[i0] = y = yn0_5
	for i1 in range(n_0 + 1, n_0 + n_1 + 1):
		k1 = h * un1_5(t, y)
		k2 = h * un1_5(t + (1/2)*h, y + (1/2)*k1)
		k3 = h * un1_5(t + (1/2)*h, y + (1/2)*k2)
		k4 = h * un1_5(t + h, y + k3)
		time1[i1] = t = time1[i1-1] + h
		RK4n_5[i1] = y = RK4n_5[i1-1] + (k1 + 2*k2 + 2*k3 + k4) / 6
	for i2 in range(n_0 + n_1 + 1, n_0 + n_1 + n_2 + 1):
		k1 = h * un2_5(t, y)
		k2 = h * un2_5(t + (1/2)*h, y + (1/2)*k1)
		k3 = h * un2_5(t + (1/2)*h, y + (1/2)*k2)
		k4 = h * un2_5(t + h, y + k3)
		time1[i2] = t = time1[i2-1] + h
		RK4n_5[i2] = y = RK4n_5[i2-1] + (k1 + 2*k2 + 2*k3 + k4) / 6
	for i3 in range(n_0 + n_1 + n_2 + 1, n_0 + n_1 + n_2 + n_3 + 1):
		k1 = h * un3_5(t, y)
		k2 = h * un3_5(t + (1/2)*h, y + (1/2)*k1)
		k3 = h * un3_5(t + (1/2)*h, y + (1/2)*k2)
		k4 = h * un3_5(t + h, y + k3)
		time1[i3] = t = time1[i3-1] + h
		RK4n_5[i3] = y = RK4n_5[i3-1] + (k1 + 2*k2 + 2*k3 + k4) / 6
	for i4 in range(n_0 + n_1 + n_2 + n_3 + 1, n_total + 1):
		k1 = h * un4_5(t, y)
		k2 = h * un4_5(t + (1/2)*h, y + (1/2)*k1)
		k3 = h * un4_5(t + (1/2)*h, y + (1/2)*k2)
		k4 = h * un4_5(t + h, y + k3)
		time1[i4] = t = time1[i4-1] + h
		RK4n_5[i4] = y = RK4n_5[i4-1] + (k1 + 2*k2 + 2*k3 + k4) / 6
	return time1, RK4n_5

# RK4 (5' released mRNA)
def rk4r_5(tt5, T1_5, T2, T3, T_last, h):
	num_total = np.arange(0, T_last, h)
	n_total = len(num_total)
	num_0 = np.arange(0, tt5, h)
	n_0 = len(num_0)
	num_1 = np.arange(tt5, T1_5, h)
	n_1 = len(num_1)
	num_2 = np.arange(T1_5, T2, h)
	n_2 = len(num_2)
	num_3 = np.arange(T2, T3, h)
	n_3 = len(num_3)
	num_4 = np.arange(T3, T_last, h)
	n_4 = len(num_4)
	time1 = [0] * (n_total + 1)
	RK4r_5 = [0] * (n_total + 1)
	time1[0] = t = 0
	RK4r_5[0] = y = yr0_5
	for i0 in range(1, n_0 + 1):
		time1[i0] = t = time1[i0-1] + h
		RK4r_5[i0] = y = yr0_5
	for i1 in range(n_0 + 1, n_0 + n_1 + 1):
		k1 = h * ur1_5(t, y)
		k2 = h * ur1_5(t + (1/2)*h, y + (1/2)*k1)
		k3 = h * ur1_5(t + (1/2)*h, y + (1/2)*k2)
		k4 = h * ur1_5(t + h, y + k3)
		time1[i1] = t = time1[i1-1] + h
		RK4r_5[i1] = y = RK4r_5[i1-1] + (k1 + 2*k2 + 2*k3 + k4) / 6
	for i2 in range(n_0 + n_1 + 1, n_0 + n_1 + n_2 + 1):
		k1 = h * ur2_5(t, y)
		k2 = h * ur2_5(t + (1/2)*h, y + (1/2)*k1)
		k3 = h * ur2_5(t + (1/2)*h, y + (1/2)*k2)
		k4 = h * ur2_5(t + h, y + k3)
		time1[i2] = t = time1[i2-1] + h
		RK4r_5[i2] = y = RK4r_5[i2-1] + (k1 + 2*k2 + 2*k3 + k4) / 6
	for i3 in range(n_0 + n_1 + n_2 + 1, n_0 + n_1 + n_2 + n_3 + 1):
		k1 = h * ur3_5(t, y)
		k2 = h * ur3_5(t + (1/2)*h, y + (1/2)*k1)
		k3 = h * ur3_5(t + (1/2)*h, y + (1/2)*k2)
		k4 = h * ur3_5(t + h, y + k3)
		time1[i3] = t = time1[i3-1] + h
		RK4r_5[i3] = y = RK4r_5[i3-1] + (k1 + 2*k2 + 2*k3 + k4) / 6
	for i4 in range(n_0 + n_1 + n_2 + n_3 + 1, n_total + 1):
		k1 = h * ur4_5(t, y)
		k2 = h * ur4_5(t + (1/2)*h, y + (1/2)*k1)
		k3 = h * ur4_5(t + (1/2)*h, y + (1/2)*k2)
		k4 = h * ur4_5(t + h, y + k3)
		time1[i4] = t = time1[i4-1] + h
		RK4r_5[i4] = y = RK4r_5[i4-1] + (k1 + 2*k2 + 2*k3 + k4) / 6
	return time1, RK4r_5

# RK4 (3' nascent mRNA)
def rk4n_3(tt3, T1_3, T2, T3, T_last, h):
	num_total = np.arange(0, T_last, h)
	n_total = len(num_total)
	num_0 = np.arange(0, tt3, h)
	n_0 = len(num_0)
	num_1 = np.arange(tt3, T2, h)
	n_1 = len(num_1)
	num_2 = np.arange(T2, T1_3, h)
	n_2 = len(num_2)
	num_3 = np.arange(T1_3, T3, h)
	n_3 = len(num_3)
	num_4 = np.arange(T3, T_last, h)
	n_4 = len(num_4)
	time1 = [0] * (n_total + 1)
	RK4n_3 = [0] * (n_total + 1)
	time1[0] = t = 0
	RK4n_3[0] = y = yn0_3
	for i0 in range(1, n_0 + 1):
		time1[i0] = t = time1[i0-1] + h
		RK4n_3[i0] = y = yn0_3
	for i1 in range(n_0 + 1, n_0 + n_1 + 1):
		k1 = h * un1_3(t, y)
		k2 = h * un1_3(t + (1/2)*h, y + (1/2)*k1)
		k3 = h * un1_3(t + (1/2)*h, y + (1/2)*k2)
		k4 = h * un1_3(t + h, y + k3)
		time1[i1] = t = time1[i1-1] + h
		RK4n_3[i1] = y = RK4n_3[i1-1] + (k1 + 2*k2 + 2*k3 + k4) / 6
	for i2 in range(n_0 + n_1 + 1, n_0 + n_1 + n_2 + 1):
		k1 = h * un2_3(t, y)
		k2 = h * un2_3(t + (1/2)*h, y + (1/2)*k1)
		k3 = h * un2_3(t + (1/2)*h, y + (1/2)*k2)
		k4 = h * un2_3(t + h, y + k3)
		time1[i2] = t = time1[i2-1] + h
		RK4n_3[i2] = y = RK4n_3[i2-1] + (k1 + 2*k2 + 2*k3 + k4) / 6
	for i3 in range(n_0 + n_1 + n_2 + 1, n_0 + n_1 + n_2 + n_3 + 1):
		k1 = h * un3_3(t, y)
		k2 = h * un3_3(t + (1/2)*h, y + (1/2)*k1)
		k3 = h * un3_3(t + (1/2)*h, y + (1/2)*k2)
		k4 = h * un3_3(t + h, y + k3)
		time1[i3] = t = time1[i3-1] + h
		RK4n_3[i3] = y = RK4n_3[i3-1] + (k1 + 2*k2 + 2*k3 + k4) / 6
	for i4 in range(n_0 + n_1 + n_2 + n_3 + 1, n_total + 1):
		k1 = h * un4_3(t, y)
		k2 = h * un4_3(t + (1/2)*h, y + (1/2)*k1)
		k3 = h * un4_3(t + (1/2)*h, y + (1/2)*k2)
		k4 = h * un4_3(t + h, y + k3)
		time1[i4] = t = time1[i4-1] + h
		RK4n_3[i4] = y = RK4n_3[i4-1] + (k1 + 2*k2 + 2*k3 + k4) / 6
	return time1, RK4n_3

# RK4 (3' released mRNA)
def rk4r_3(tt3, T1_3, T2, T3, T_last, h):
	num_total = np.arange(0, T_last, h)
	n_total = len(num_total)
	num_0 = np.arange(0, tt3, h)
	n_0 = len(num_0)
	num_1 = np.arange(tt3, T2, h)
	n_1 = len(num_1)
	num_2 = np.arange(T2, T1_3, h)
	n_2 = len(num_2)
	num_3 = np.arange(T1_3, T3, h)
	n_3 = len(num_3)
	num_4 = np.arange(T3, T_last, h)
	n_4 = len(num_4)
	time1 = [0] * (n_total + 1)
	RK4r_3 = [0] * (n_total + 1)
	time1[0] = t = 0
	RK4r_3[0] = y = yr0_3
	for i0 in range(1, n_0 + 1):
		time1[i0] = t = time1[i0-1] + h
		RK4r_3[i0] = y = yr0_3
	for i1 in range(n_0 + 1, n_0 + n_1 + 1):
		k1 = h * ur1_3(t, y)
		k2 = h * ur1_3(t + (1/2)*h, y + (1/2)*k1)
		k3 = h * ur1_3(t + (1/2)*h, y + (1/2)*k2)
		k4 = h * ur1_3(t + h, y + k3)
		time1[i1] = t = time1[i1-1] + h
		RK4r_3[i1] = y = RK4r_3[i1-1] + (k1 + 2*k2 + 2*k3 + k4) / 6
	for i2 in range(n_0 + n_1 + 1, n_0 + n_1 + n_2 + 1):
		k1 = h * ur2_3(t, y)
		k2 = h * ur2_3(t + (1/2)*h, y + (1/2)*k1)
		k3 = h * ur2_3(t + (1/2)*h, y + (1/2)*k2)
		k4 = h * ur2_3(t + h, y + k3)
		time1[i2] = t = time1[i2-1] + h
		RK4r_3[i2] = y = RK4r_3[i2-1] + (k1 + 2*k2 + 2*k3 + k4) / 6
	for i3 in range(n_0 + n_1 + n_2 + 1, n_0 + n_1 + n_2 + n_3 + 1):
		k1 = h * ur3_3(t, y)
		k2 = h * ur3_3(t + (1/2)*h, y + (1/2)*k1)
		k3 = h * ur3_3(t + (1/2)*h, y + (1/2)*k2)
		k4 = h * ur3_3(t + h, y + k3)
		time1[i3] = t = time1[i3-1] + h
		RK4r_3[i3] = y = RK4r_3[i3-1] + (k1 + 2*k2 + 2*k3 + k4) / 6
	for i4 in range(n_0 + n_1 + n_2 + n_3 + 1, n_total + 1):
		k1 = h * ur4_3(t, y)
		k2 = h * ur4_3(t + (1/2)*h, y + (1/2)*k1)
		k3 = h * ur4_3(t + (1/2)*h, y + (1/2)*k2)
		k4 = h * ur4_3(t + h, y + k3)
		time1[i4] = t = time1[i4-1] + h
		RK4r_3[i4] = y = RK4r_3[i4-1] + (k1 + 2*k2 + 2*k3 + k4) / 6
	return time1, RK4r_3

# Calculation
time1, RK4n_5 = rk4n_5(tt5, T1_5, T2, T3, 600, 0.1)
time1, RK4r_5 = rk4r_5(tt5, T1_5, T2, T3, 600, 0.1)
time1, RK4n_3 = rk4n_3(tt3, T1_3, T2, T3, 600, 0.1)
time1, RK4r_3 = rk4r_3(tt3, T1_3, T2, T3, 600, 0.1)

def RK4_5():
	n_total = len(np.arange(0, 600, 0.1))
	RK4sum_5 = [0] * (n_total + 1)
	for i in range (0, n_total + 1):
		RK4sum_5[i] = RK4n_5[i] + RK4r_5[i]
	return RK4sum_5

def RK4_3():
	n_total = len(np.arange(0, 600, 0.1))
	RK4sum_3 = [0] * (n_total + 1)
	for i in range (0, n_total + 1):
		RK4sum_3[i] = RK4n_3[i] + RK4r_3[i]
	return RK4sum_3

# Calculation
RK4sum_5 = RK4_5()
RK4sum_3 = RK4_3()

# Columns in output excel file
col1 = "time [s]"
col2 = "RK4_n (5')"
col3 = "RK4_r (5')"
col4 = "RK4_n (3')"
col5 = "RK4_r (3')"
col6 = "RK4 n+r (5')"
col7 = "RK4 n+r (3')"

# Output excel file
data = DataFrame({col1:time1, col2:RK4n_5, col3:RK4r_5, col4:RK4n_3, col5:RK4r_3, col6:RK4sum_5, col7:RK4sum_3})
data.to_excel('case1.xlsx', sheet_name='sheet1', index=False)
