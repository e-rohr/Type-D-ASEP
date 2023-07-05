import re

eDual5122 = "r**(-6) * (-q**5*F1*F2*F3*F4*F3*F5*F2 + q**4*F1*F2*F3*F5*F2*F3*F4 - q**3*F1*F2*F3*F5*F4*F2*F3 + q**6*F1*F2*F3*F5*F4*F3*F2 + q**4*F1*F2*F4*F3*F5*F2*F3 - q**5*F1*F2*F5*F3*F4*F3*F2 - (q**4 - q**2)*F1*F2**2*F3*F5*F4*F3 + q**4*F1*F3*F2*F4*F3*F5*F2 - q**5*F1*F3*F2*F5*F4*F3*F2 + q**4*F1*F5*F3*F2*F4*F3*F2 + (q**3 - q)*F2*F1*F2*F3*F5*F4*F3 - q**3*F2*F1*F3*F4*F2*F3*F5 + (q**4 + q**2)*F2*F1*F3*F4*F3*F5*F2 - q**3*F2*F1*F3*F5*F2*F3*F4 + (q**4 + q**2)*F2*F1*F3*F5*F4*F2*F3 - (q**5 + q)*F2*F1*F3*F5*F4*F3*F2 - q**3*F2*F1*F4*F3*F5*F2*F3 - q**3*F2*F1*F5*F3*F4*F2*F3 + (q**4 + q**2)*F2*F1*F5*F3*F4*F3*F2 - q*F2*F3*F4*F3*F5*F2*F1 + q**2*F2*F3*F5*F2*F1*F3*F4 - q**3*F2*F3*F5*F4*F2*F1*F3 + 1*F2*F3*F5*F4*F3*F2*F1 + q**2*F2*F4*F3*F5*F2*F1*F3 - q*F2*F5*F3*F4*F3*F2*F1 - q**3*F3*F2*F1*F4*F3*F5*F2 - q**3*F3*F2*F1*F5*F3*F4*F2 + (q**4 + q**2)*F3*F2*F1*F5*F4*F3*F2 + q**2*F3*F2*F4*F3*F5*F2*F1 - q*F3*F2*F5*F4*F3*F2*F1 - q**3*F4*F3*F2*F1*F5*F3*F2 - q**3*F5*F3*F2*F1*F4*F3*F2 + q**2*F5*F3*F2*F4*F3*F2*F1)"
eDual5132 = "r**(-5) * (-q**4*F1*F2*F3*F4*F3*F5 + q**5*F1*F2*F3*F5*F4*F3 - q**4*F1*F2*F5*F3*F4*F3 + q**3*F1*F3*F2*F4*F3*F5 + q**3*F1*F3*F2*F5*F3*F4 - (q**4 + q**2)*F1*F3*F2*F5*F4*F3 - q**2*F1*F3*F4*F3*F2*F5 + q*F1*F3*F5*F4*F3*F2 + q**3*F1*F4*F3*F2*F5*F3 + q**3*F1*F5*F3*F2*F4*F3 - q**2*F1*F5*F3*F4*F3*F2 + q**3*F2*F1*F3*F4*F3*F5 - q**4*F2*F1*F3*F5*F4*F3 + q**3*F2*F1*F5*F3*F4*F3 - q**2*F3*F2*F1*F4*F3*F5 - q**2*F3*F2*F1*F5*F3*F4 + (q**3 + q)*F3*F2*F1*F5*F4*F3 + q*F3*F4*F3*F2*F1*F5 - 1*F3*F5*F4*F3*F2*F1 - q**2*F4*F3*F2*F1*F5*F3 - q**2*F5*F3*F2*F1*F4*F3 + q*F5*F3*F4*F3*F2*F1)"

test = "F1*F2*F3*F4*F3*F5*F2"


def reverse(match_obj):
    l = len(match_obj.groups())
    str = ""
    for i in range(1, l + 1):
        str += "F"
        str += match_obj.group(l + 1 - i)
        if (i < l):
            str += "*"
    return str

eDual5312 = re.sub(r"F(\d)\*F(\d)\*F(\d)\*F(\d)\*F(\d)\*F(\d)", reverse, eDual5132)

print(eDual5312)

