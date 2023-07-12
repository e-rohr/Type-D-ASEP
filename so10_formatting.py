S = "\frac{r^{2} F_{4} K_{4} KNeg_{5} KNeg_{4}^{2} K_{4} E_{4}}{q} + \frac{r^{4} \left(\frac{q F_{3} F_{4}}{r} - \frac{F_{4} F_{3}}{r}\right) K_{4} KNeg_{5} KNeg_{3}^{2} KNeg_{4}^{2} K_{3} K_{4} \left(\frac{q E_{4} E_{3}}{r} - \frac{E_{3} E_{4}}{r}\right)}{q^{3}} + \frac{r^{2} F_{3} K_{4} KNeg_{5} KNeg_{3}^{2} KNeg_{4}^{2} K_{3} E_{3}}{q^{3}} + \frac{r^{6} \left(\frac{q^{2} F_{2} F_{3} F_{4}}{r^{2}} - \frac{q F_{2} F_{4} F_{3}}{r^{2}} - \frac{q F_{3} F_{2} F_{4}}{r^{2}} + \frac{F_{4} F_{3} F_{2}}{r^{2}}\right) K_{4} KNeg_{5} KNeg_{2}^{2} KNeg_{3}^{2} KNeg_{4}^{2} K_{2} K_{3} K_{4} \left(\frac{q^{2} E_{4} E_{3} E_{2}}{r^{2}} - \frac{q E_{3} E_{4} E_{2}}{r^{2}} - \frac{q E_{4} E_{2} E_{3}}{r^{2}} + \frac{E_{2} E_{3} E_{4}}{r^{2}}\right)}{q^{5}} + \frac{r^{4} \left(\frac{q F_{2} F_{3}}{r} - \frac{F_{3} F_{2}}{r}\right) K_{4} KNeg_{5} KNeg_{2}^{2} KNeg_{3}^{2} KNeg_{4}^{2} K_{2} K_{3} \left(\frac{q E_{3} E_{2}}{r} - \frac{E_{2} E_{3}}{r}\right)}{q^{5}} + \frac{r^{2} F_{2} K_{4} KNeg_{5} KNeg_{2}^{2} KNeg_{3}^{2} KNeg_{4}^{2} K_{2} E_{2}}{q^{5}} + \frac{r^{8} \left(\frac{q^{3} F_{1} F_{2} F_{3} F_{4}}{r^{3}} - \frac{q^{2} F_{1} F_{2} F_{4} F_{3}}{r^{3}} - \frac{q^{2} F_{1} F_{3} F_{2} F_{4}}{r^{3}} - \frac{q^{2} F_{2} F_{1} F_{3} F_{4}}{r^{3}} + \frac{q F_{1} F_{4} F_{3} F_{2}}{r^{3}} + \frac{q F_{2} F_{1} F_{4} F_{3}}{r^{3}} + \frac{q F_{3} F_{2} F_{1} F_{4}}{r^{3}} - \frac{F_{4} F_{3} F_{2} F_{1}}{r^{3}}\right) K_{4} KNeg_{5} KNeg_{1}^{2} KNeg_{2}^{2} KNeg_{3}^{2} KNeg_{4}^{2} K_{1} K_{2} K_{3} K_{4} \left(\frac{q^{3} E_{4} E_{3} E_{2} E_{1}}{r^{3}} - \frac{q^{2} E_{3} E_{4} E_{2} E_{1}}{r^{3}} - \frac{q^{2} E_{4} E_{2} E_{3} E_{1}}{r^{3}} - \frac{q^{2} E_{4} E_{3} E_{1} E_{2}}{r^{3}} + \frac{q E_{2} E_{3} E_{4} E_{1}}{r^{3}} + \frac{q E_{3} E_{4} E_{1} E_{2}}{r^{3}} + \frac{q E_{4} E_{1} E_{2} E_{3}}{r^{3}} - \frac{E_{1} E_{2} E_{3} E_{4}}{r^{3}}\right)}{q^{7}} + \frac{r^{6} \left(\frac{q^{2} F_{1} F_{2} F_{3}}{r^{2}} - \frac{q F_{1} F_{3} F_{2}}{r^{2}} - \frac{q F_{2} F_{1} F_{3}}{r^{2}} + \frac{F_{3} F_{2} F_{1}}{r^{2}}\right) K_{4} KNeg_{5} KNeg_{1}^{2} KNeg_{2}^{2} KNeg_{3}^{2} KNeg_{4}^{2} K_{1} K_{2} K_{3} \left(\frac{q^{2} E_{3} E_{2} E_{1}}{r^{2}} - \frac{q E_{2} E_{3} E_{1}}{r^{2}} - \frac{q E_{3} E_{1} E_{2}}{r^{2}} + \frac{E_{1} E_{2} E_{3}}{r^{2}}\right)}{q^{7}} + \frac{r^{4} \left(\frac{q F_{1} F_{2}}{r} - \frac{F_{2} F_{1}}{r}\right) K_{4} KNeg_{5} KNeg_{1}^{2} KNeg_{2}^{2} KNeg_{3}^{2} KNeg_{4}^{2} K_{1} K_{2} \left(\frac{q E_{2} E_{1}}{r} - \frac{E_{1} E_{2}}{r}\right)}{q^{7}} + \frac{r^{2} F_{1} K_{4} KNeg_{5} KNeg_{1}^{2} KNeg_{2}^{2} KNeg_{3}^{2} KNeg_{4}^{2} K_{1} E_{1}}{q^{7}}"

for i in range(len(S)):
    if S[i]=="E" or S[i] == "F":
        start = i
        indices = ""
        count = 0
        while S[i] == "E" or S[i] == "F":
            i += 3
            count += 2
            indices += S[i]
            i += 3
            count += 2
        end = i
        S = S[:i+2] + indices + S[end - 1:]

print(S)