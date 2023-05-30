print("hi")

if member.BESS.P_max > 0
    P_min = -P_max
if P_t1 > 0
    alpha = P_t1
else
    alpha = P_t1 - P_min
end

if P_t2 >= 0
    beta = P_max - P_t2
else
    beta = -P_t2 
end

if t1 <= t2
    e = minimum(member.SOC[t1:t2])
else
    e = member.SOC_max - maximum(member.SOC[t1:t2])
end