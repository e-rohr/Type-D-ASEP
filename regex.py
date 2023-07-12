import re
#import regex




input = r"""
- q r^{8} \left(\frac{q^{3} F_{5} F_{3} F_{2} F_{1}}{r^{3}} - \frac{q^{2} F_{3} F_{5} F_{2} F_{1}}{r^{3}} - \frac{q^{2} F_{5} F_{2} F_{3} F_{1}}{r^{3}} - \frac{q^{2} F_{5} F_{3} F_{1} F_{2}}{r^{3}} + \frac{q F_{2} F_{3} F_{5} F_{1}}{r^{3}} + \frac{q F_{3} F_{5} F_{1} F_{2}}{r^{3}} + \frac{q F_{5} F_{1} F_{2} F_{3}}{r^{3}} - \frac{F_{1} F_{2} F_{3} F_{5}}{r^{3}}\right) K_{1} K_{2} K_{3} K_{4} \left(\frac{q^{3} E_{1} E_{2} E_{3} E_{5}}{r^{3}} - \frac{q^{2} E_{1} E_{2} E_{5} E_{3}}{r^{3}} - \frac{q^{2} E_{1} E_{3} E_{2} E_{5}}{r^{3}} - \frac{q^{2} E_{2} E_{1} E_{3} E_{5}}{r^{3}} + \frac{q E_{1} E_{5} E_{3} E_{2}}{r^{3}} + \frac{q E_{2} E_{1} E_{5} E_{3}}{r^{3}} + \frac{q E_{3} E_{2} E_{1} E_{5}}{r^{3}} - \frac{E_{5} E_{3} E_{2} E_{1}}{r^{3}}\right) - q r^{6} \left(\frac{q^{2} F_{5} F_{3} F_{2}}{r^{2}} - \frac{q F_{3} F_{5} F_{2}}{r^{2}} - \frac{q F_{5} F_{2} F_{3}}{r^{2}} + \frac{F_{2} F_{3} F_{5}}{r^{2}}\right) K_{2} K_{3} K_{4} \left(\frac{q^{2} E_{2} E_{3} E_{5}}{r^{2}} - \frac{q E_{2} E_{5} E_{3}}{r^{2}} - \frac{q E_{3} E_{2} E_{5}}{r^{2}} + \frac{E_{5} E_{3} E_{2}}{r^{2}}\right) - q r^{4} \left(\frac{q F_{5} F_{3}}{r} - \frac{F_{3} F_{5}}{r}\right) K_{3} K_{4} \left(\frac{q E_{3} E_{5}}{r} - \frac{E_{5} E_{3}}{r}\right) - q r^{2} F_{5} K_{4} E_{5} - r^{4} F_{5} F_{4} E_{5} E_{4} - \frac{r^{10} \left(\frac{q^{4} F_{4} F_{5} F_{3} F_{2} F_{1}}{r^{4}} - \frac{q^{3} F_{4} F_{3} F_{5} F_{2} F_{1}}{r^{4}} - \frac{q^{3} F_{4} F_{5} F_{2} F_{3} F_{1}}{r^{4}} - \frac{q^{3} F_{4} F_{5} F_{3} F_{1} F_{2}}{r^{4}} - \frac{q^{3} F_{5} F_{3} F_{4} F_{2} F_{1}}{r^{4}} + \frac{q^{2} F_{3} F_{4} F_{5} F_{2} F_{1}}{r^{4}} + \frac{q^{2} F_{4} F_{2} F_{3} F_{5} F_{1}}{r^{4}} + \frac{q^{2} F_{4} F_{3} F_{5} F_{1} F_{2}}{r^{4}} + \frac{q^{2} F_{4} F_{5} F_{1} F_{2} F_{3}}{r^{4}} + \frac{q^{2} F_{5} F_{2} F_{3} F_{4} F_{1}}{r^{4}} + \frac{q^{2} F_{5} F_{3} F_{4} F_{1} F_{2}}{r^{4}} - \frac{q F_{2} F_{3} F_{4} F_{5} F_{1}}{r^{4}} - \frac{q F_{3} F_{4} F_{5} F_{1} F_{2}}{r^{4}} - \frac{q F_{4} F_{1} F_{2} F_{3} F_{5}}{r^{4}} - \frac{q F_{5} F_{1} F_{2} F_{3} F_{4}}{r^{4}} + \frac{F_{1} F_{2} F_{3} F_{4} F_{5}}{r^{4}}\right) K_{1} K_{2} K_{3} \left(\frac{q^{4} E_{1} E_{2} E_{3} E_{5} E_{4}}{r^{4}} - \frac{q^{3} E_{1} E_{2} E_{4} E_{3} E_{5}}{r^{4}} - \frac{q^{3} E_{1} E_{2} E_{5} E_{3} E_{4}}{r^{4}} - \frac{q^{3} E_{1} E_{3} E_{2} E_{5} E_{4}}{r^{4}} - \frac{q^{3} E_{2} E_{1} E_{3} E_{5} E_{4}}{r^{4}} + \frac{q^{2} E_{1} E_{2} E_{5} E_{4} E_{3}}{r^{4}} + \frac{q^{2} E_{1} E_{4} E_{3} E_{2} E_{5}}{r^{4}} + \frac{q^{2} E_{1} E_{5} E_{3} E_{2} E_{4}}{r^{4}} + \frac{q^{2} E_{2} E_{1} E_{4} E_{3} E_{5}}{r^{4}} + \frac{q^{2} E_{2} E_{1} E_{5} E_{3} E_{4}}{r^{4}} + \frac{q^{2} E_{3} E_{2} E_{1} E_{5} E_{4}}{r^{4}} - \frac{q E_{1} E_{5} E_{4} E_{3} E_{2}}{r^{4}} - \frac{q E_{2} E_{1} E_{5} E_{4} E_{3}}{r^{4}} - \frac{q E_{4} E_{3} E_{2} E_{1} E_{5}}{r^{4}} - \frac{q E_{5} E_{3} E_{2} E_{1} E_{4}}{r^{4}} + \frac{E_{5} E_{4} E_{3} E_{2} E_{1}}{r^{4}}\right)}{q} - \frac{r^{8} \left(\frac{q^{3} F_{4} F_{5} F_{3} F_{2}}{r^{3}} - \frac{q^{2} F_{4} F_{3} F_{5} F_{2}}{r^{3}} - \frac{q^{2} F_{4} F_{5} F_{2} F_{3}}{r^{3}} - \frac{q^{2} F_{5} F_{3} F_{4} F_{2}}{r^{3}} + \frac{q F_{3} F_{4} F_{5} F_{2}}{r^{3}} + \frac{q F_{4} F_{2} F_{3} F_{5}}{r^{3}} + \frac{q F_{5} F_{2} F_{3} F_{4}}{r^{3}} - \frac{F_{2} F_{3} F_{4} F_{5}}{r^{3}}\right) K_{2} K_{3} \left(\frac{q^{3} E_{2} E_{3} E_{5} E_{4}}{r^{3}} - \frac{q^{2} E_{2} E_{4} E_{3} E_{5}}{r^{3}} - \frac{q^{2} E_{2} E_{5} E_{3} E_{4}}{r^{3}} - \frac{q^{2} E_{3} E_{2} E_{5} E_{4}}{r^{3}} + \frac{q E_{2} E_{5} E_{4} E_{3}}{r^{3}} + \frac{q E_{4} E_{3} E_{2} E_{5}}{r^{3}} + \frac{q E_{5} E_{3} E_{2} E_{4}}{r^{3}} - \frac{E_{5} E_{4} E_{3} E_{2}}{r^{3}}\right)}{q} - \frac{r^{6} \left(\frac{q^{2} F_{4} F_{5} F_{3}}{r^{2}} - \frac{q F_{4} F_{3} F_{5}}{r^{2}} - \frac{q F_{5} F_{3} F_{4}}{r^{2}} + \frac{F_{3} F_{4} F_{5}}{r^{2}}\right) K_{3} \left(\frac{q^{2} E_{3} E_{5} E_{4}}{r^{2}} - \frac{q E_{4} E_{3} E_{5}}{r^{2}} - \frac{q E_{5} E_{3} E_{4}}{r^{2}} + \frac{E_{5} E_{4} E_{3}}{r^{2}}\right)}{q} - \frac{r^{2} F_{5} KNeg_{4} E_{5}}{q} - \frac{r^{8} \left(- \frac{q F_{3} F_{4} F_{3} F_{5}}{r^{2}} - \frac{q F_{5} F_{3} F_{4} F_{3}}{r^{2}} + \frac{\left(q^{4} + q^{2}\right) F_{3} F_{5} F_{4} F_{3}}{q^{2} r^{2}}\right) \left(- \frac{q E_{3} E_{4} E_{3} E_{5}}{r^{2}} - \frac{q E_{5} E_{3} E_{4} E_{3}}{r^{2}} + \frac{\left(q^{4} + q^{2}\right) E_{3} E_{5} E_{4} E_{3}}{q^{2} r^{2}}\right)}{q^{2}} - \frac{r^{10} \left(\frac{q^{4} F_{3} F_{4} F_{5} F_{3} F_{2}}{r^{4}} - \frac{q^{2} F_{3} F_{5} F_{3} F_{4} F_{2}}{r^{3}} - \frac{q^{2} F_{4} F_{3} F_{5} F_{3} F_{2}}{r^{3}} + \frac{q^{2} F_{3} F_{4} F_{2} F_{3} F_{5}}{r^{4}} + \frac{q^{2} F_{3} F_{5} F_{2} F_{3} F_{4}}{r^{4}} + \frac{q^{2} F_{4} F_{3} F_{5} F_{2} F_{3}}{r^{4}} + \frac{q^{2} F_{5} F_{3} F_{4} F_{2} F_{3}}{r^{4}} - \frac{q F_{3} F_{2} F_{3} F_{4} F_{5}}{r^{4}} - \frac{q F_{4} F_{3} F_{2} F_{3} F_{5}}{r^{4}} - \frac{q F_{4} F_{5} F_{3} F_{2} F_{3}}{r^{4}} - \frac{q F_{5} F_{3} F_{2} F_{3} F_{4}}{r^{4}} + \frac{F_{2} F_{3} F_{4} F_{5} F_{3}}{r^{4}} + \frac{\left(- q^{7} - q^{5}\right) F_{3} F_{4} F_{5} F_{2} F_{3}}{q^{4} r^{4}}\right) K_{2} \left(\frac{q^{4} E_{2} E_{3} E_{5} E_{4} E_{3}}{r^{4}} - \frac{q^{3} E_{2} E_{3} E_{4} E_{3} E_{5}}{r^{4}} - \frac{q^{3} E_{2} E_{5} E_{3} E_{4} E_{3}}{r^{4}} + \frac{q^{2} E_{3} E_{2} E_{4} E_{3} E_{5}}{r^{4}} + \frac{q^{2} E_{3} E_{2} E_{5} E_{3} E_{4}}{r^{4}} + \frac{q^{2} E_{4} E_{3} E_{2} E_{5} E_{3}}{r^{4}} + \frac{q^{2} E_{5} E_{3} E_{2} E_{4} E_{3}}{r^{4}} - \frac{q E_{3} E_{4} E_{3} E_{2} E_{5}}{r^{4}} - \frac{q E_{5} E_{3} E_{4} E_{3} E_{2}}{r^{4}} + \frac{E_{3} E_{5} E_{4} E_{3} E_{2}}{r^{4}} + \frac{\left(- q^{7} - q^{5}\right) E_{3} E_{2} E_{5} E_{4} E_{3}}{q^{4} r^{4}}\right)}{q^{3}} - \frac{r^{6} \left(\frac{q^{2} F_{3} F_{5} F_{4}}{r^{2}} - \frac{q F_{4} F_{3} F_{5}}{r^{2}} - \frac{q F_{5} F_{3} F_{4}}{r^{2}} + \frac{F_{5} F_{4} F_{3}}{r^{2}}\right) KNeg_{3} \left(\frac{q^{2} E_{4} E_{5} E_{3}}{r^{2}} - \frac{q E_{4} E_{3} E_{5}}{r^{2}} - \frac{q E_{5} E_{3} E_{4}}{r^{2}} + \frac{E_{3} E_{4} E_{5}}{r^{2}}\right)}{q^{3}} - \frac{r^{4} \left(\frac{q F_{3} F_{5}}{r} - \frac{F_{5} F_{3}}{r}\right) KNeg_{3} KNeg_{4} \left(\frac{q E_{5} E_{3}}{r} - \frac{E_{3} E_{5}}{r}\right)}{q^{3}} - \frac{r^{2} \left(q^{5} F_{3} F_{5} F_{4} F_{3} F_{2} F_{1} - q^{4} F_{3} F_{5} F_{4} F_{3} F_{1} F_{2} - q^{3} r F_{3} F_{4} F_{3} F_{5} F_{2} F_{1} - q^{3} r F_{5} F_{3} F_{4} F_{3} F_{2} F_{1} + q^{3} F_{3} F_{4} F_{2} F_{3} F_{5} F_{1} + q^{3} F_{3} F_{5} F_{2} F_{3} F_{4} F_{1} + q^{3} F_{4} F_{3} F_{5} F_{2} F_{3} F_{1} + q^{3} F_{5} F_{3} F_{4} F_{2} F_{3} F_{1} - q^{2} F_{3} F_{2} F_{3} F_{5} F_{4} F_{1} - q^{2} F_{3} F_{4} F_{1} F_{2} F_{3} F_{5} - q^{2} F_{3} F_{5} F_{1} F_{2} F_{3} F_{4} - q^{2} F_{4} F_{3} F_{2} F_{3} F_{5} F_{1} - q^{2} F_{4} F_{3} F_{5} F_{1} F_{2} F_{3} - q^{2} F_{5} F_{3} F_{2} F_{3} F_{4} F_{1} - q^{2} F_{5} F_{3} F_{4} F_{1} F_{2} F_{3} - q^{2} F_{5} F_{4} F_{3} F_{2} F_{3} F_{1} + q F_{2} F_{3} F_{5} F_{4} F_{3} F_{1} + q F_{3} F_{1} F_{2} F_{3} F_{5} F_{4} + q F_{4} F_{3} F_{1} F_{2} F_{3} F_{5} + q F_{5} F_{3} F_{1} F_{2} F_{3} F_{4} + q F_{5} F_{4} F_{3} F_{1} F_{2} F_{3} + \left(q^{3} - q\right) F_{3} F_{4} F_{3} F_{5} F_{1} F_{2} + \left(q^{3} - q\right) F_{5} F_{3} F_{4} F_{3} F_{1} F_{2} + \left(q^{3} + q\right) F_{3} F_{5} F_{4} F_{1} F_{2} F_{3} + \left(- q^{4} - q^{2}\right) F_{3} F_{5} F_{4} F_{2} F_{3} F_{1} - F_{1} F_{2} F_{3} F_{5} F_{4} F_{3}\right) K_{1} K_{2} \left(q^{5} E_{1} E_{2} E_{3} E_{5} E_{4} E_{3} - q^{4} E_{1} E_{2} E_{3} E_{4} E_{3} E_{5} - q^{4} E_{1} E_{2} E_{5} E_{3} E_{4} E_{3} - q^{4} E_{2} E_{1} E_{3} E_{5} E_{4} E_{3} + q^{3} E_{1} E_{3} E_{2} E_{4} E_{3} E_{5} + q^{3} E_{1} E_{3} E_{2} E_{5} E_{3} E_{4} + q^{3} E_{1} E_{4} E_{3} E_{2} E_{5} E_{3} + q^{3} E_{1} E_{5} E_{3} E_{2} E_{4} E_{3} + q^{3} E_{2} E_{1} E_{3} E_{4} E_{3} E_{5} + q^{3} E_{2} E_{1} E_{5} E_{3} E_{4} E_{3} - q^{2} E_{1} E_{3} E_{4} E_{3} E_{2} E_{5} - q^{2} E_{1} E_{5} E_{3} E_{4} E_{3} E_{2} - q^{2} E_{3} E_{2} E_{1} E_{4} E_{3} E_{5} - q^{2} E_{3} E_{2} E_{1} E_{5} E_{3} E_{4} - q^{2} E_{4} E_{3} E_{2} E_{1} E_{5} E_{3} - q^{2} E_{5} E_{3} E_{2} E_{1} E_{4} E_{3} + q E_{1} E_{3} E_{5} E_{4} E_{3} E_{2} + q E_{3} E_{4} E_{3} E_{2} E_{1} E_{5} + q E_{5} E_{3} E_{4} E_{3} E_{2} E_{1} + \left(q^{3} + q\right) E_{3} E_{2} E_{1} E_{5} E_{4} E_{3} + \left(- q^{4} - q^{2}\right) E_{1} E_{3} E_{2} E_{5} E_{4} E_{3} - E_{3} E_{5} E_{4} E_{3} E_{2} E_{1}\right)}{q^{3}} - \frac{r^{12} \left(\frac{q^{2} F_{2} F_{3} F_{5} F_{2} F_{3} F_{4}}{r^{4}} + \frac{q^{2} F_{2} F_{4} F_{3} F_{5} F_{2} F_{3}}{r^{4}} - \frac{q^{2} F_{2}^{2} F_{3} F_{5} F_{4} F_{3}}{r^{4}} + \frac{q^{2} F_{3} F_{2} F_{4} F_{3} F_{5} F_{2}}{r^{4}} + \frac{q^{2} F_{5} F_{3} F_{2} F_{4} F_{3} F_{2}}{r^{4}} + \frac{\left(- q^{7} - q^{5}\right) F_{2} F_{3} F_{4} F_{3} F_{5} F_{2}}{q^{4} r^{4}} + \frac{\left(- q^{7} - q^{5}\right) F_{2} F_{5} F_{3} F_{4} F_{3} F_{2}}{q^{4} r^{4}} + \frac{\left(- q^{7} - q^{5}\right) F_{3} F_{2} F_{5} F_{4} F_{3} F_{2}}{q^{4} r^{4}} + \frac{\left(q^{8} + q^{6} + q^{4}\right) F_{2} F_{3} F_{5} F_{4} F_{3} F_{2}}{q^{4} r^{4}}\right) \left(\frac{q^{2} E_{2} E_{3} E_{5} E_{2} E_{3} E_{4}}{r^{4}} + \frac{q^{2} E_{2} E_{4} E_{3} E_{5} E_{2} E_{3}}{r^{4}} - \frac{q^{2} E_{2}^{2} E_{3} E_{5} E_{4} E_{3}}{r^{4}} + \frac{q^{2} E_{3} E_{2} E_{4} E_{3} E_{5} E_{2}}{r^{4}} + \frac{q^{2} E_{5} E_{3} E_{2} E_{4} E_{3} E_{2}}{r^{4}} + \frac{\left(- q^{7} - q^{5}\right) E_{2} E_{3} E_{4} E_{3} E_{5} E_{2}}{q^{4} r^{4}} + \frac{\left(- q^{7} - q^{5}\right) E_{2} E_{5} E_{3} E_{4} E_{3} E_{2}}{q^{4} r^{4}} + \frac{\left(- q^{7} - q^{5}\right) E_{3} E_{2} E_{5} E_{4} E_{3} E_{2}}{q^{4} r^{4}} + \frac{\left(q^{8} + q^{6} + q^{4}\right) E_{2} E_{3} E_{5} E_{4} E_{3} E_{2}}{q^{4} r^{4}}\right)}{q^{4}} - \frac{r^{10} \left(\frac{q^{4} F_{2} F_{3} F_{5} F_{4} F_{3}}{r^{4}} - \frac{q^{3} F_{2} F_{3} F_{4} F_{3} F_{5}}{r^{4}} - \frac{q^{3} F_{2} F_{5} F_{3} F_{4} F_{3}}{r^{4}} + \frac{q^{2} F_{3} F_{2} F_{4} F_{3} F_{5}}{r^{4}} + \frac{q^{2} F_{3} F_{2} F_{5} F_{3} F_{4}}{r^{4}} + \frac{q^{2} F_{4} F_{3} F_{2} F_{5} F_{3}}{r^{4}} + \frac{q^{2} F_{5} F_{3} F_{2} F_{4} F_{3}}{r^{4}} - \frac{q F_{3} F_{4} F_{3} F_{2} F_{5}}{r^{4}} - \frac{q F_{5} F_{3} F_{4} F_{3} F_{2}}{r^{4}} + \frac{F_{3} F_{5} F_{4} F_{3} F_{2}}{r^{4}} + \frac{\left(- q^{7} - q^{5}\right) F_{3} F_{2} F_{5} F_{4} F_{3}}{q^{4} r^{4}}\right) KNeg_{2} \left(\frac{q^{4} E_{3} E_{4} E_{5} E_{3} E_{2}}{r^{4}} - \frac{q^{2} E_{3} E_{5} E_{3} E_{4} E_{2}}{r^{3}} - \frac{q^{2} E_{4} E_{3} E_{5} E_{3} E_{2}}{r^{3}} + \frac{q^{2} E_{3} E_{4} E_{2} E_{3} E_{5}}{r^{4}} + \frac{q^{2} E_{3} E_{5} E_{2} E_{3} E_{4}}{r^{4}} + \frac{q^{2} E_{4} E_{3} E_{5} E_{2} E_{3}}{r^{4}} + \frac{q^{2} E_{5} E_{3} E_{4} E_{2} E_{3}}{r^{4}} - \frac{q E_{3} E_{2} E_{3} E_{4} E_{5}}{r^{4}} - \frac{q E_{4} E_{3} E_{2} E_{3} E_{5}}{r^{4}} - \frac{q E_{4} E_{5} E_{3} E_{2} E_{3}}{r^{4}} - \frac{q E_{5} E_{3} E_{2} E_{3} E_{4}}{r^{4}} + \frac{E_{2} E_{3} E_{4} E_{5} E_{3}}{r^{4}} + \frac{\left(- q^{7} - q^{5}\right) E_{3} E_{4} E_{5} E_{2} E_{3}}{q^{4} r^{4}}\right)}{q^{5}} - \frac{r^{8} \left(\frac{q^{3} F_{2} F_{3} F_{5} F_{4}}{r^{3}} - \frac{q^{2} F_{2} F_{4} F_{3} F_{5}}{r^{3}} - \frac{q^{2} F_{2} F_{5} F_{3} F_{4}}{r^{3}} - \frac{q^{2} F_{3} F_{2} F_{5} F_{4}}{r^{3}} + \frac{q F_{2} F_{5} F_{4} F_{3}}{r^{3}} + \frac{q F_{4} F_{3} F_{2} F_{5}}{r^{3}} + \frac{q F_{5} F_{3} F_{2} F_{4}}{r^{3}} - \frac{F_{5} F_{4} F_{3} F_{2}}{r^{3}}\right) KNeg_{2} KNeg_{3} \left(\frac{q^{3} E_{4} E_{5} E_{3} E_{2}}{r^{3}} - \frac{q^{2} E_{4} E_{3} E_{5} E_{2}}{r^{3}} - \frac{q^{2} E_{4} E_{5} E_{2} E_{3}}{r^{3}} - \frac{q^{2} E_{5} E_{3} E_{4} E_{2}}{r^{3}} + \frac{q E_{3} E_{4} E_{5} E_{2}}{r^{3}} + \frac{q E_{4} E_{2} E_{3} E_{5}}{r^{3}} + \frac{q E_{5} E_{2} E_{3} E_{4}}{r^{3}} - \frac{E_{2} E_{3} E_{4} E_{5}}{r^{3}}\right)}{q^{5}} - \frac{r^{6} \left(\frac{q^{2} F_{2} F_{3} F_{5}}{r^{2}} - \frac{q F_{2} F_{5} F_{3}}{r^{2}} - \frac{q F_{3} F_{2} F_{5}}{r^{2}} + \frac{F_{5} F_{3} F_{2}}{r^{2}}\right) KNeg_{2} KNeg_{3} KNeg_{4} \left(\frac{q^{2} E_{5} E_{3} E_{2}}{r^{2}} - \frac{q E_{3} E_{5} E_{2}}{r^{2}} - \frac{q E_{5} E_{2} E_{3}}{r^{2}} + \frac{E_{2} E_{3} E_{5}}{r^{2}}\right)}{q^{5}} - \frac{r^{2} \left(q^{6} F_{2} F_{3} F_{5} F_{4} F_{3} F_{2} F_{1} + q^{4} F_{2} F_{3} F_{4} F_{3} F_{5} F_{1} F_{2} + q^{4} F_{2} F_{5} F_{3} F_{4} F_{3} F_{1} F_{2} + q^{3} r F_{2} F_{3} F_{5} F_{2} F_{3} F_{4} F_{1} + q^{3} r F_{2} F_{4} F_{3} F_{5} F_{2} F_{3} F_{1} - q^{3} r F_{2}^{2} F_{3} F_{5} F_{4} F_{3} F_{1} + q^{3} r F_{3} F_{2} F_{4} F_{3} F_{5} F_{2} F_{1} + q^{3} r F_{5} F_{3} F_{2} F_{4} F_{3} F_{2} F_{1} - q^{3} F_{2} F_{3} F_{4} F_{1} F_{2} F_{3} F_{5} - q^{3} F_{2} F_{3} F_{5} F_{1} F_{2} F_{3} F_{4} - q^{3} F_{2} F_{4} F_{3} F_{5} F_{1} F_{2} F_{3} - q^{3} F_{2} F_{5} F_{3} F_{4} F_{1} F_{2} F_{3} - q^{3} F_{3} F_{2} F_{4} F_{3} F_{5} F_{1} F_{2} - q^{3} F_{3} F_{2} F_{5} F_{3} F_{4} F_{1} F_{2} - q^{3} F_{4} F_{3} F_{2} F_{5} F_{3} F_{1} F_{2} - q^{3} F_{5} F_{3} F_{2} F_{4} F_{3} F_{1} F_{2} + q^{2} F_{2} F_{3} F_{1} F_{2} F_{3} F_{5} F_{4} + q^{2} F_{2} F_{4} F_{3} F_{1} F_{2} F_{3} F_{5} + q^{2} F_{2} F_{5} F_{3} F_{1} F_{2} F_{3} F_{4} + q^{2} F_{2} F_{5} F_{4} F_{3} F_{1} F_{2} F_{3} + q^{2} F_{3} F_{2} F_{4} F_{1} F_{2} F_{3} F_{5} + q^{2} F_{3} F_{2} F_{5} F_{1} F_{2} F_{3} F_{4} + q^{2} F_{3} F_{4} F_{3} F_{2} F_{5} F_{1} F_{2} + q^{2} F_{4} F_{3} F_{2} F_{5} F_{1} F_{2} F_{3} + q^{2} F_{5} F_{3} F_{2} F_{4} F_{1} F_{2} F_{3} + q^{2} F_{5} F_{3} F_{4} F_{3} F_{2} F_{1} F_{2} - q F_{2} F_{1} F_{2} F_{3} F_{5} F_{4} F_{3} - q F_{3} F_{2} F_{1} F_{2} F_{3} F_{5} F_{4} - q F_{3} F_{5} F_{4} F_{3} F_{2} F_{1} F_{2} - q F_{4} F_{3} F_{2} F_{1} F_{2} F_{3} F_{5} - q F_{5} F_{3} F_{2} F_{1} F_{2} F_{3} F_{4} - q F_{5} F_{4} F_{3} F_{2} F_{1} F_{2} F_{3} + \left(- q^{3} - q\right) F_{3} F_{2} F_{5} F_{4} F_{1} F_{2} F_{3} + \left(q^{4} + q^{2}\right) F_{2} F_{3} F_{5} F_{4} F_{1} F_{2} F_{3} + \left(q^{4} + q^{2}\right) F_{3} F_{2} F_{5} F_{4} F_{3} F_{1} F_{2} + \left(- q^{5} - q\right) F_{2} F_{3} F_{5} F_{4} F_{3} F_{1} F_{2} + \left(- q^{5} + q\right) F_{2} F_{3} F_{4} F_{3} F_{5} F_{2} F_{1} + \left(- q^{5} + q\right) F_{2} F_{5} F_{3} F_{4} F_{3} F_{2} F_{1} + \left(- q^{5} + q\right) F_{3} F_{2} F_{5} F_{4} F_{3} F_{2} F_{1} + F_{1} F_{2} F_{3} F_{5} F_{4} F_{3} F_{2}\right) K_{1} \left(q^{6} E_{1} E_{2} E_{3} E_{5} E_{4} E_{3} E_{2} - q^{5} E_{1} E_{2} E_{3} E_{4} E_{3} E_{5} E_{2} - q^{5} E_{1} E_{2} E_{5} E_{3} E_{4} E_{3} E_{2} - q^{5} E_{1} E_{3} E_{2} E_{5} E_{4} E_{3} E_{2} + q^{4} E_{1} E_{2} E_{3} E_{5} E_{2} E_{3} E_{4} + q^{4} E_{1} E_{2} E_{4} E_{3} E_{5} E_{2} E_{3} + q^{4} E_{1} E_{3} E_{2} E_{4} E_{3} E_{5} E_{2} + q^{4} E_{1} E_{5} E_{3} E_{2} E_{4} E_{3} E_{2} - q^{3} r E_{1} E_{2}^{2} E_{3} E_{5} E_{4} E_{3} - q^{3} E_{1} E_{2} E_{3} E_{5} E_{4} E_{2} E_{3} - q^{3} E_{2} E_{1} E_{3} E_{4} E_{2} E_{3} E_{5} - q^{3} E_{2} E_{1} E_{3} E_{5} E_{2} E_{3} E_{4} - q^{3} E_{2} E_{1} E_{4} E_{3} E_{5} E_{2} E_{3} - q^{3} E_{2} E_{1} E_{5} E_{3} E_{4} E_{2} E_{3} - q^{3} E_{2} E_{3} E_{5} E_{4} E_{2} E_{1} E_{3} - q^{3} E_{3} E_{2} E_{1} E_{4} E_{3} E_{5} E_{2} - q^{3} E_{3} E_{2} E_{1} E_{5} E_{3} E_{4} E_{2} - q^{3} E_{4} E_{3} E_{2} E_{1} E_{5} E_{3} E_{2} - q^{3} E_{5} E_{3} E_{2} E_{1} E_{4} E_{3} E_{2} + q^{2} E_{2} E_{3} E_{5} E_{2} E_{1} E_{3} E_{4} + q^{2} E_{2} E_{4} E_{3} E_{5} E_{2} E_{1} E_{3} + q^{2} E_{3} E_{2} E_{4} E_{3} E_{5} E_{2} E_{1} + q^{2} E_{5} E_{3} E_{2} E_{4} E_{3} E_{2} E_{1} - q E_{2} E_{3} E_{4} E_{3} E_{5} E_{2} E_{1} - q E_{2} E_{5} E_{3} E_{4} E_{3} E_{2} E_{1} - q E_{3} E_{2} E_{5} E_{4} E_{3} E_{2} E_{1} + \left(q^{3} - q\right) E_{2} E_{1} E_{2} E_{3} E_{5} E_{4} E_{3} + \left(q^{4} + q^{2}\right) E_{2} E_{1} E_{3} E_{4} E_{3} E_{5} E_{2} + \left(q^{4} + q^{2}\right) E_{2} E_{1} E_{3} E_{5} E_{4} E_{2} E_{3} + \left(q^{4} + q^{2}\right) E_{2} E_{1} E_{5} E_{3} E_{4} E_{3} E_{2} + \left(q^{4} + q^{2}\right) E_{3} E_{2} E_{1} E_{5} E_{4} E_{3} E_{2} + \left(- q^{5} - q\right) E_{2} E_{1} E_{3} E_{5} E_{4} E_{3} E_{2} + E_{2} E_{3} E_{5} E_{4} E_{3} E_{2} E_{1}\right)}{q^{5}} - \frac{r^{4} \left(- q^{3} F_{1} F_{2} F_{3} F_{4} F_{1} F_{2} F_{3} F_{5} + q^{3} F_{1} F_{2} F_{3} F_{5} F_{4} F_{3} F_{1} F_{2} - q^{3} F_{1} F_{2} F_{5} F_{3} F_{4} F_{1} F_{2} F_{3} + q^{3} F_{1} F_{3} F_{2} F_{1} F_{2} F_{3} F_{5} F_{4} - q^{3} F_{1} F_{3} F_{2} F_{5} F_{3} F_{4} F_{1} F_{2} + q^{3} F_{1} F_{3} F_{5} F_{4} F_{3} F_{2} F_{1} F_{2} + q^{3} F_{1} F_{4} F_{3} F_{2} F_{1} F_{2} F_{3} F_{5} - q^{3} F_{1} F_{4} F_{3} F_{2} F_{5} F_{3} F_{1} F_{2} + q^{3} F_{1} F_{5} F_{3} F_{2} F_{1} F_{2} F_{3} F_{4} + q^{3} F_{1} F_{5} F_{4} F_{3} F_{2} F_{1} F_{2} F_{3} + q^{3} F_{2} F_{1} F_{2} F_{3} F_{5} F_{4} F_{3} F_{1} - q^{3} F_{2} F_{1} F_{3} F_{5} F_{2} F_{3} F_{4} F_{1} - q^{3} F_{2} F_{1} F_{4} F_{3} F_{5} F_{2} F_{3} F_{1} - q^{3} F_{3} F_{2} F_{1} F_{4} F_{3} F_{5} F_{2} F_{1} - q^{3} F_{5} F_{3} F_{2} F_{1} F_{4} F_{3} F_{2} F_{1} + q^{2} \left(- q r - 2\right) F_{1} F_{2}^{2} F_{3} F_{5} F_{4} F_{3} F_{1} + q^{2} \left(- q r - 2\right) F_{1}^{2} F_{2} F_{3} F_{5} F_{4} F_{3} F_{2} + q^{2} \left(q r + 2\right) F_{1} F_{2} F_{3} F_{5} F_{2} F_{3} F_{4} F_{1} + q^{2} \left(q r + 2\right) F_{1} F_{2} F_{3} F_{5} F_{4} F_{1} F_{2} F_{3} + q^{2} \left(q r + 2\right) F_{1} F_{2} F_{4} F_{3} F_{5} F_{2} F_{3} F_{1} + q^{2} \left(q r + 2\right) F_{1} F_{3} F_{2} F_{4} F_{3} F_{5} F_{2} F_{1} + q^{2} \left(q r + 2\right) F_{1} F_{5} F_{3} F_{2} F_{4} F_{3} F_{2} F_{1} + q^{2} \left(q r + 2\right) F_{2} F_{1} F_{3} F_{4} F_{3} F_{5} F_{2} F_{1} + q^{2} \left(q r + 2\right) F_{2} F_{1} F_{5} F_{3} F_{4} F_{3} F_{2} F_{1} + q^{2} \left(q r + 2\right) F_{3} F_{2} F_{1} F_{5} F_{4} F_{3} F_{2} F_{1} - q \left(q r + 2\right)^{2} F_{1} F_{2} F_{3} F_{4} F_{3} F_{5} F_{2} F_{1} - q \left(q r + 2\right)^{2} F_{1} F_{2} F_{5} F_{3} F_{4} F_{3} F_{2} F_{1} - q \left(q r + 2\right)^{2} F_{1} F_{3} F_{2} F_{5} F_{4} F_{3} F_{2} F_{1} - q \left(q^{4} + q^{2} + 1\right) F_{2} F_{1} F_{3} F_{5} F_{4} F_{3} F_{2} F_{1} + \left(q^{6} + q^{4} + q^{2} + 1\right) F_{1} F_{2} F_{3} F_{5} F_{4} F_{3} F_{2} F_{1}\right) \left(- q^{3} E_{1} E_{2} E_{3} E_{4} E_{1} E_{2} E_{3} E_{5} + q^{3} E_{1} E_{2} E_{3} E_{5} E_{4} E_{3} E_{1} E_{2} - q^{3} E_{1} E_{2} E_{5} E_{3} E_{4} E_{1} E_{2} E_{3} + q^{3} E_{1} E_{3} E_{2} E_{1} E_{2} E_{3} E_{5} E_{4} - q^{3} E_{1} E_{3} E_{2} E_{5} E_{3} E_{4} E_{1} E_{2} + q^{3} E_{1} E_{3} E_{5} E_{4} E_{3} E_{2} E_{1} E_{2} + q^{3} E_{1} E_{4} E_{3} E_{2} E_{1} E_{2} E_{3} E_{5} - q^{3} E_{1} E_{4} E_{3} E_{2} E_{5} E_{3} E_{1} E_{2} + q^{3} E_{1} E_{5} E_{3} E_{2} E_{1} E_{2} E_{3} E_{4} + q^{3} E_{1} E_{5} E_{4} E_{3} E_{2} E_{1} E_{2} E_{3} + q^{3} E_{2} E_{1} E_{2} E_{3} E_{5} E_{4} E_{3} E_{1} - q^{3} E_{2} E_{1} E_{3} E_{5} E_{2} E_{3} E_{4} E_{1} - q^{3} E_{2} E_{1} E_{4} E_{3} E_{5} E_{2} E_{3} E_{1} - q^{3} E_{3} E_{2} E_{1} E_{4} E_{3} E_{5} E_{2} E_{1} - q^{3} E_{5} E_{3} E_{2} E_{1} E_{4} E_{3} E_{2} E_{1} + q^{2} \left(- q r - 2\right) E_{1} E_{2}^{2} E_{3} E_{5} E_{4} E_{3} E_{1} + q^{2} \left(- q r - 2\right) E_{1}^{2} E_{2} E_{3} E_{5} E_{4} E_{3} E_{2} + q^{2} \left(q r + 2\right) E_{1} E_{2} E_{3} E_{5} E_{2} E_{3} E_{4} E_{1} + q^{2} \left(q r + 2\right) E_{1} E_{2} E_{3} E_{5} E_{4} E_{1} E_{2} E_{3} + q^{2} \left(q r + 2\right) E_{1} E_{2} E_{4} E_{3} E_{5} E_{2} E_{3} E_{1} + q^{2} \left(q r + 2\right) E_{1} E_{3} E_{2} E_{4} E_{3} E_{5} E_{2} E_{1} + q^{2} \left(q r + 2\right) E_{1} E_{5} E_{3} E_{2} E_{4} E_{3} E_{2} E_{1} + q^{2} \left(q r + 2\right) E_{2} E_{1} E_{3} E_{4} E_{3} E_{5} E_{2} E_{1} + q^{2} \left(q r + 2\right) E_{2} E_{1} E_{5} E_{3} E_{4} E_{3} E_{2} E_{1} + q^{2} \left(q r + 2\right) E_{3} E_{2} E_{1} E_{5} E_{4} E_{3} E_{2} E_{1} - q \left(q r + 2\right)^{2} E_{1} E_{2} E_{3} E_{4} E_{3} E_{5} E_{2} E_{1} - q \left(q r + 2\right)^{2} E_{1} E_{2} E_{5} E_{3} E_{4} E_{3} E_{2} E_{1} - q \left(q r + 2\right)^{2} E_{1} E_{3} E_{2} E_{5} E_{4} E_{3} E_{2} E_{1} - q \left(q^{4} + q^{2} + 1\right) E_{2} E_{1} E_{3} E_{5} E_{4} E_{3} E_{2} E_{1} + \left(q^{6} + q^{4} + q^{2} + 1\right) E_{1} E_{2} E_{3} E_{5} E_{4} E_{3} E_{2} E_{1}\right)}{q^{6}} - \frac{r^{10} \left(\frac{q^{4} F_{1} F_{2} F_{3} F_{5} F_{4}}{r^{4}} - \frac{q^{3} F_{1} F_{2} F_{4} F_{3} F_{5}}{r^{4}} - \frac{q^{3} F_{1} F_{2} F_{5} F_{3} F_{4}}{r^{4}} - \frac{q^{3} F_{1} F_{3} F_{2} F_{5} F_{4}}{r^{4}} - \frac{q^{3} F_{2} F_{1} F_{3} F_{5} F_{4}}{r^{4}} + \frac{q^{2} F_{1} F_{2} F_{5} F_{4} F_{3}}{r^{4}} + \frac{q^{2} F_{1} F_{4} F_{3} F_{2} F_{5}}{r^{4}} + \frac{q^{2} F_{1} F_{5} F_{3} F_{2} F_{4}}{r^{4}} + \frac{q^{2} F_{2} F_{1} F_{4} F_{3} F_{5}}{r^{4}} + \frac{q^{2} F_{2} F_{1} F_{5} F_{3} F_{4}}{r^{4}} + \frac{q^{2} F_{3} F_{2} F_{1} F_{5} F_{4}}{r^{4}} - \frac{q F_{1} F_{5} F_{4} F_{3} F_{2}}{r^{4}} - \frac{q F_{2} F_{1} F_{5} F_{4} F_{3}}{r^{4}} - \frac{q F_{4} F_{3} F_{2} F_{1} F_{5}}{r^{4}} - \frac{q F_{5} F_{3} F_{2} F_{1} F_{4}}{r^{4}} + \frac{F_{5} F_{4} F_{3} F_{2} F_{1}}{r^{4}}\right) KNeg_{1} KNeg_{2} KNeg_{3} \left(\frac{q^{4} E_{4} E_{5} E_{3} E_{2} E_{1}}{r^{4}} - \frac{q^{3} E_{4} E_{3} E_{5} E_{2} E_{1}}{r^{4}} - \frac{q^{3} E_{4} E_{5} E_{2} E_{3} E_{1}}{r^{4}} - \frac{q^{3} E_{4} E_{5} E_{3} E_{1} E_{2}}{r^{4}} - \frac{q^{3} E_{5} E_{3} E_{4} E_{2} E_{1}}{r^{4}} + \frac{q^{2} E_{3} E_{4} E_{5} E_{2} E_{1}}{r^{4}} + \frac{q^{2} E_{4} E_{2} E_{3} E_{5} E_{1}}{r^{4}} + \frac{q^{2} E_{4} E_{3} E_{5} E_{1} E_{2}}{r^{4}} + \frac{q^{2} E_{4} E_{5} E_{1} E_{2} E_{3}}{r^{4}} + \frac{q^{2} E_{5} E_{2} E_{3} E_{4} E_{1}}{r^{4}} + \frac{q^{2} E_{5} E_{3} E_{4} E_{1} E_{2}}{r^{4}} - \frac{q E_{2} E_{3} E_{4} E_{5} E_{1}}{r^{4}} - \frac{q E_{3} E_{4} E_{5} E_{1} E_{2}}{r^{4}} - \frac{q E_{4} E_{1} E_{2} E_{3} E_{5}}{r^{4}} - \frac{q E_{5} E_{1} E_{2} E_{3} E_{4}}{r^{4}} + \frac{E_{1} E_{2} E_{3} E_{4} E_{5}}{r^{4}}\right)}{q^{7}} - \frac{r^{8} \left(\frac{q^{3} F_{1} F_{2} F_{3} F_{5}}{r^{3}} - \frac{q^{2} F_{1} F_{2} F_{5} F_{3}}{r^{3}} - \frac{q^{2} F_{1} F_{3} F_{2} F_{5}}{r^{3}} - \frac{q^{2} F_{2} F_{1} F_{3} F_{5}}{r^{3}} + \frac{q F_{1} F_{5} F_{3} F_{2}}{r^{3}} + \frac{q F_{2} F_{1} F_{5} F_{3}}{r^{3}} + \frac{q F_{3} F_{2} F_{1} F_{5}}{r^{3}} - \frac{F_{5} F_{3} F_{2} F_{1}}{r^{3}}\right) KNeg_{1} KNeg_{2} KNeg_{3} KNeg_{4} \left(\frac{q^{3} E_{5} E_{3} E_{2} E_{1}}{r^{3}} - \frac{q^{2} E_{3} E_{5} E_{2} E_{1}}{r^{3}} - \frac{q^{2} E_{5} E_{2} E_{3} E_{1}}{r^{3}} - \frac{q^{2} E_{5} E_{3} E_{1} E_{2}}{r^{3}} + \frac{q E_{2} E_{3} E_{5} E_{1}}{r^{3}} + \frac{q E_{3} E_{5} E_{1} E_{2}}{r^{3}} + \frac{q E_{5} E_{1} E_{2} E_{3}}{r^{3}} - \frac{E_{1} E_{2} E_{3} E_{5}}{r^{3}}\right)}{q^{7}} - \frac{r^{2} \left(q^{5} F_{1} F_{2} F_{3} F_{5} F_{4} F_{3} - q^{4} F_{1} F_{2} F_{3} F_{4} F_{3} F_{5} - q^{4} F_{1} F_{2} F_{5} F_{3} F_{4} F_{3} - q^{4} F_{2} F_{1} F_{3} F_{5} F_{4} F_{3} + q^{3} F_{1} F_{3} F_{2} F_{4} F_{3} F_{5} + q^{3} F_{1} F_{3} F_{2} F_{5} F_{3} F_{4} + q^{3} F_{1} F_{4} F_{3} F_{2} F_{5} F_{3} + q^{3} F_{1} F_{5} F_{3} F_{2} F_{4} F_{3} + q^{3} F_{2} F_{1} F_{3} F_{4} F_{3} F_{5} + q^{3} F_{2} F_{1} F_{5} F_{3} F_{4} F_{3} - q^{2} F_{1} F_{3} F_{4} F_{3} F_{2} F_{5} - q^{2} F_{1} F_{5} F_{3} F_{4} F_{3} F_{2} - q^{2} F_{3} F_{2} F_{1} F_{4} F_{3} F_{5} - q^{2} F_{3} F_{2} F_{1} F_{5} F_{3} F_{4} - q^{2} F_{4} F_{3} F_{2} F_{1} F_{5} F_{3} - q^{2} F_{5} F_{3} F_{2} F_{1} F_{4} F_{3} + q F_{1} F_{3} F_{5} F_{4} F_{3} F_{2} + q F_{3} F_{4} F_{3} F_{2} F_{1} F_{5} + q F_{5} F_{3} F_{4} F_{3} F_{2} F_{1} + \left(q^{3} + q\right) F_{3} F_{2} F_{1} F_{5} F_{4} F_{3} + \left(- q^{4} - q^{2}\right) F_{1} F_{3} F_{2} F_{5} F_{4} F_{3} - F_{3} F_{5} F_{4} F_{3} F_{2} F_{1}\right) KNeg_{1} KNeg_{2} \left(q^{5} E_{3} E_{5} E_{4} E_{3} E_{2} E_{1} - q^{4} E_{3} E_{5} E_{4} E_{3} E_{1} E_{2} - q^{3} r E_{3} E_{4} E_{3} E_{5} E_{2} E_{1} - q^{3} r E_{5} E_{3} E_{4} E_{3} E_{2} E_{1} + q^{3} E_{3} E_{4} E_{2} E_{3} E_{5} E_{1} + q^{3} E_{3} E_{5} E_{2} E_{3} E_{4} E_{1} + q^{3} E_{4} E_{3} E_{5} E_{2} E_{3} E_{1} + q^{3} E_{5} E_{3} E_{4} E_{2} E_{3} E_{1} - q^{2} E_{3} E_{2} E_{3} E_{5} E_{4} E_{1} - q^{2} E_{3} E_{4} E_{1} E_{2} E_{3} E_{5} - q^{2} E_{3} E_{5} E_{1} E_{2} E_{3} E_{4} - q^{2} E_{4} E_{3} E_{2} E_{3} E_{5} E_{1} - q^{2} E_{4} E_{3} E_{5} E_{1} E_{2} E_{3} - q^{2} E_{5} E_{3} E_{2} E_{3} E_{4} E_{1} - q^{2} E_{5} E_{3} E_{4} E_{1} E_{2} E_{3} - q^{2} E_{5} E_{4} E_{3} E_{2} E_{3} E_{1} + q E_{2} E_{3} E_{5} E_{4} E_{3} E_{1} + q E_{3} E_{1} E_{2} E_{3} E_{5} E_{4} + q E_{4} E_{3} E_{1} E_{2} E_{3} E_{5} + q E_{5} E_{3} E_{1} E_{2} E_{3} E_{4} + q E_{5} E_{4} E_{3} E_{1} E_{2} E_{3} + \left(q^{3} - q\right) E_{3} E_{4} E_{3} E_{5} E_{1} E_{2} + \left(q^{3} - q\right) E_{5} E_{3} E_{4} E_{3} E_{1} E_{2} + \left(q^{3} + q\right) E_{3} E_{5} E_{4} E_{1} E_{2} E_{3} + \left(- q^{4} - q^{2}\right) E_{3} E_{5} E_{4} E_{2} E_{3} E_{1} - E_{1} E_{2} E_{3} E_{5} E_{4} E_{3}\right)}{q^{7}} - \frac{r^{2} \left(q^{6} F_{1} F_{2} F_{3} F_{5} F_{4} F_{3} F_{2} - q^{5} F_{1} F_{2} F_{3} F_{4} F_{3} F_{5} F_{2} - q^{5} F_{1} F_{2} F_{5} F_{3} F_{4} F_{3} F_{2} - q^{5} F_{1} F_{3} F_{2} F_{5} F_{4} F_{3} F_{2} + q^{4} F_{1} F_{2} F_{3} F_{5} F_{2} F_{3} F_{4} + q^{4} F_{1} F_{2} F_{4} F_{3} F_{5} F_{2} F_{3} + q^{4} F_{1} F_{3} F_{2} F_{4} F_{3} F_{5} F_{2} + q^{4} F_{1} F_{5} F_{3} F_{2} F_{4} F_{3} F_{2} - q^{3} r F_{1} F_{2}^{2} F_{3} F_{5} F_{4} F_{3} - q^{3} F_{1} F_{2} F_{3} F_{5} F_{4} F_{2} F_{3} - q^{3} F_{2} F_{1} F_{3} F_{4} F_{2} F_{3} F_{5} - q^{3} F_{2} F_{1} F_{3} F_{5} F_{2} F_{3} F_{4} - q^{3} F_{2} F_{1} F_{4} F_{3} F_{5} F_{2} F_{3} - q^{3} F_{2} F_{1} F_{5} F_{3} F_{4} F_{2} F_{3} - q^{3} F_{2} F_{3} F_{5} F_{4} F_{2} F_{1} F_{3} - q^{3} F_{3} F_{2} F_{1} F_{4} F_{3} F_{5} F_{2} - q^{3} F_{3} F_{2} F_{1} F_{5} F_{3} F_{4} F_{2} - q^{3} F_{4} F_{3} F_{2} F_{1} F_{5} F_{3} F_{2} - q^{3} F_{5} F_{3} F_{2} F_{1} F_{4} F_{3} F_{2} + q^{2} F_{2} F_{3} F_{5} F_{2} F_{1} F_{3} F_{4} + q^{2} F_{2} F_{4} F_{3} F_{5} F_{2} F_{1} F_{3} + q^{2} F_{3} F_{2} F_{4} F_{3} F_{5} F_{2} F_{1} + q^{2} F_{5} F_{3} F_{2} F_{4} F_{3} F_{2} F_{1} - q F_{2} F_{3} F_{4} F_{3} F_{5} F_{2} F_{1} - q F_{2} F_{5} F_{3} F_{4} F_{3} F_{2} F_{1} - q F_{3} F_{2} F_{5} F_{4} F_{3} F_{2} F_{1} + \left(q^{3} - q\right) F_{2} F_{1} F_{2} F_{3} F_{5} F_{4} F_{3} + \left(q^{4} + q^{2}\right) F_{2} F_{1} F_{3} F_{4} F_{3} F_{5} F_{2} + \left(q^{4} + q^{2}\right) F_{2} F_{1} F_{3} F_{5} F_{4} F_{2} F_{3} + \left(q^{4} + q^{2}\right) F_{2} F_{1} F_{5} F_{3} F_{4} F_{3} F_{2} + \left(q^{4} + q^{2}\right) F_{3} F_{2} F_{1} F_{5} F_{4} F_{3} F_{2} + \left(- q^{5} - q\right) F_{2} F_{1} F_{3} F_{5} F_{4} F_{3} F_{2} + F_{2} F_{3} F_{5} F_{4} F_{3} F_{2} F_{1}\right) KNeg_{1} \left(q^{6} E_{2} E_{3} E_{5} E_{4} E_{3} E_{2} E_{1} + q^{4} E_{2} E_{3} E_{4} E_{3} E_{5} E_{1} E_{2} + q^{4} E_{2} E_{5} E_{3} E_{4} E_{3} E_{1} E_{2} + q^{3} r E_{2} E_{3} E_{5} E_{2} E_{3} E_{4} E_{1} + q^{3} r E_{2} E_{4} E_{3} E_{5} E_{2} E_{3} E_{1} - q^{3} r E_{2}^{2} E_{3} E_{5} E_{4} E_{3} E_{1} + q^{3} r E_{3} E_{2} E_{4} E_{3} E_{5} E_{2} E_{1} + q^{3} r E_{5} E_{3} E_{2} E_{4} E_{3} E_{2} E_{1} - q^{3} E_{2} E_{3} E_{4} E_{1} E_{2} E_{3} E_{5} - q^{3} E_{2} E_{3} E_{5} E_{1} E_{2} E_{3} E_{4} - q^{3} E_{2} E_{4} E_{3} E_{5} E_{1} E_{2} E_{3} - q^{3} E_{2} E_{5} E_{3} E_{4} E_{1} E_{2} E_{3} - q^{3} E_{3} E_{2} E_{4} E_{3} E_{5} E_{1} E_{2} - q^{3} E_{3} E_{2} E_{5} E_{3} E_{4} E_{1} E_{2} - q^{3} E_{4} E_{3} E_{2} E_{5} E_{3} E_{1} E_{2} - q^{3} E_{5} E_{3} E_{2} E_{4} E_{3} E_{1} E_{2} + q^{2} E_{2} E_{3} E_{1} E_{2} E_{3} E_{5} E_{4} + q^{2} E_{2} E_{4} E_{3} E_{1} E_{2} E_{3} E_{5} + q^{2} E_{2} E_{5} E_{3} E_{1} E_{2} E_{3} E_{4} + q^{2} E_{2} E_{5} E_{4} E_{3} E_{1} E_{2} E_{3} + q^{2} E_{3} E_{2} E_{4} E_{1} E_{2} E_{3} E_{5} + q^{2} E_{3} E_{2} E_{5} E_{1} E_{2} E_{3} E_{4} + q^{2} E_{3} E_{4} E_{3} E_{2} E_{5} E_{1} E_{2} + q^{2} E_{4} E_{3} E_{2} E_{5} E_{1} E_{2} E_{3} + q^{2} E_{5} E_{3} E_{2} E_{4} E_{1} E_{2} E_{3} + q^{2} E_{5} E_{3} E_{4} E_{3} E_{2} E_{1} E_{2} - q E_{2} E_{1} E_{2} E_{3} E_{5} E_{4} E_{3} - q E_{3} E_{2} E_{1} E_{2} E_{3} E_{5} E_{4} - q E_{3} E_{5} E_{4} E_{3} E_{2} E_{1} E_{2} - q E_{4} E_{3} E_{2} E_{1} E_{2} E_{3} E_{5} - q E_{5} E_{3} E_{2} E_{1} E_{2} E_{3} E_{4} - q E_{5} E_{4} E_{3} E_{2} E_{1} E_{2} E_{3} + \left(- q^{3} - q\right) E_{3} E_{2} E_{5} E_{4} E_{1} E_{2} E_{3} + \left(q^{4} + q^{2}\right) E_{2} E_{3} E_{5} E_{4} E_{1} E_{2} E_{3} + \left(q^{4} + q^{2}\right) E_{3} E_{2} E_{5} E_{4} E_{3} E_{1} E_{2} + \left(- q^{5} - q\right) E_{2} E_{3} E_{5} E_{4} E_{3} E_{1} E_{2} + \left(- q^{5} + q\right) E_{2} E_{3} E_{4} E_{3} E_{5} E_{2} E_{1} + \left(- q^{5} + q\right) E_{2} E_{5} E_{3} E_{4} E_{3} E_{2} E_{1} + \left(- q^{5} + q\right) E_{3} E_{2} E_{5} E_{4} E_{3} E_{2} E_{1} + E_{1} E_{2} E_{3} E_{5} E_{4} E_{3} E_{2}\right)}{q^{7}}
"""


#eDual5122 = "r**(-6) * (-q**5*F1*F2*F3*F4*F3*F5*F2 + q**4*F1*F2*F3*F5*F2*F3*F4 - q**3*F1*F2*F3*F5*F4*F2*F3 + q**6*F1*F2*F3*F5*F4*F3*F2 + q**4*F1*F2*F4*F3*F5*F2*F3 - q**5*F1*F2*F5*F3*F4*F3*F2 - (q**4 - q**2)*F1*F2**2*F3*F5*F4*F3 + q**4*F1*F3*F2*F4*F3*F5*F2 - q**5*F1*F3*F2*F5*F4*F3*F2 + q**4*F1*F5*F3*F2*F4*F3*F2 + (q**3 - q)*F2*F1*F2*F3*F5*F4*F3 - q**3*F2*F1*F3*F4*F2*F3*F5 + (q**4 + q**2)*F2*F1*F3*F4*F3*F5*F2 - q**3*F2*F1*F3*F5*F2*F3*F4 + (q**4 + q**2)*F2*F1*F3*F5*F4*F2*F3 - (q**5 + q)*F2*F1*F3*F5*F4*F3*F2 - q**3*F2*F1*F4*F3*F5*F2*F3 - q**3*F2*F1*F5*F3*F4*F2*F3 + (q**4 + q**2)*F2*F1*F5*F3*F4*F3*F2 - q*F2*F3*F4*F3*F5*F2*F1 + q**2*F2*F3*F5*F2*F1*F3*F4 - q**3*F2*F3*F5*F4*F2*F1*F3 + 1*F2*F3*F5*F4*F3*F2*F1 + q**2*F2*F4*F3*F5*F2*F1*F3 - q*F2*F5*F3*F4*F3*F2*F1 - q**3*F3*F2*F1*F4*F3*F5*F2 - q**3*F3*F2*F1*F5*F3*F4*F2 + (q**4 + q**2)*F3*F2*F1*F5*F4*F3*F2 + q**2*F3*F2*F4*F3*F5*F2*F1 - q*F3*F2*F5*F4*F3*F2*F1 - q**3*F4*F3*F2*F1*F5*F3*F2 - q**3*F5*F3*F2*F1*F4*F3*F2 + q**2*F5*F3*F2*F4*F3*F2*F1)"
#eDual5132 = "r**(-5) * (-q**4*F1*F2*F3*F4*F3*F5 + q**5*F1*F2*F3*F5*F4*F3 - q**4*F1*F2*F5*F3*F4*F3 + q**3*F1*F3*F2*F4*F3*F5 + q**3*F1*F3*F2*F5*F3*F4 - (q**4 + q**2)*F1*F3*F2*F5*F4*F3 - q**2*F1*F3*F4*F3*F2*F5 + q*F1*F3*F5*F4*F3*F2 + q**3*F1*F4*F3*F2*F5*F3 + q**3*F1*F5*F3*F2*F4*F3 - q**2*F1*F5*F3*F4*F3*F2 + q**3*F2*F1*F3*F4*F3*F5 - q**4*F2*F1*F3*F5*F4*F3 + q**3*F2*F1*F5*F3*F4*F3 - q**2*F3*F2*F1*F4*F3*F5 - q**2*F3*F2*F1*F5*F3*F4 + (q**3 + q)*F3*F2*F1*F5*F4*F3 + q*F3*F4*F3*F2*F1*F5 - 1*F3*F5*F4*F3*F2*F1 - q**2*F4*F3*F2*F1*F5*F3 - q**2*F5*F3*F2*F1*F4*F3 + q*F5*F3*F4*F3*F2*F1)"
#eDual5112 = "r**(-6) * (-q**3*F1*F2*F3*F4*F1*F2*F3*F5 - (q**5 + 2*q**3 + q)*F1*F2*F3*F4*F3*F5*F2*F1 + (q**4 + q**2)*F1*F2*F3*F5*F2*F3*F4*F1 + (q**4 + q**2)*F1*F2*F3*F5*F4*F1*F2*F3 + q**3*F1*F2*F3*F5*F4*F3*F1*F2 + (q**6 + q**4 + q**2 + 1)*F1*F2*F3*F5*F4*F3*F2*F1 + (q**4 + q**2)*F1*F2*F4*F3*F5*F2*F3*F1 - q**3*F1*F2*F5*F3*F4*F1*F2*F3 - (q**5 + 2*q**3 + q)*F1*F2*F5*F3*F4*F3*F2*F1 - (q**4 + q**2)*F1*F2**2*F3*F5*F4*F3*F1 + q**3*F1*F3*F2*F1*F2*F3*F5*F4 + (q**4 + q**2)*F1*F3*F2*F4*F3*F5*F2*F1 - q**3*F1*F3*F2*F5*F3*F4*F1*F2 - (q**5 + 2*q**3 + q)*F1*F3*F2*F5*F4*F3*F2*F1 + q**3*F1*F3*F5*F4*F3*F2*F1*F2 + q**3*F1*F4*F3*F2*F1*F2*F3*F5 - q**3*F1*F4*F3*F2*F5*F3*F1*F2 + q**3*F1*F5*F3*F2*F1*F2*F3*F4 + (q**4 + q**2)*F1*F5*F3*F2*F4*F3*F2*F1 + q**3*F1*F5*F4*F3*F2*F1*F2*F3 - (q**4 + q**2)*F1**2*F2*F3*F5*F4*F3*F2 + q**3*F2*F1*F2*F3*F5*F4*F3*F1 + (q**4 + q**2)*F2*F1*F3*F4*F3*F5*F2*F1 - q**3*F2*F1*F3*F5*F2*F3*F4*F1 - (q**5 + q**3 + q)*F2*F1*F3*F5*F4*F3*F2*F1 - q**3*F2*F1*F4*F3*F5*F2*F3*F1 + (q**4 + q**2)*F2*F1*F5*F3*F4*F3*F2*F1 - q**3*F3*F2*F1*F4*F3*F5*F2*F1 + (q**4 + q**2)*F3*F2*F1*F5*F4*F3*F2*F1 - q**3*F5*F3*F2*F1*F4*F3*F2*F1)"

#test = "F_{1}*F_{2}*F_{3}*F_{4}*F_{5}*F_{6}*F_{7}*F_{8}"




def combineF(match_obj):
    l = len(match_obj.groups())
    str = "F_{"
    for i in range(1, l+1):
        str += match_obj.group(i)
    return str + "}"

def combineE(match_obj):
    l = len(match_obj.groups())
    str = "E_{"
    for i in range(1, l+1):
        str += match_obj.group(i)
    return str + "}"

def combineK(match_obj):
    l = len(match_obj.groups())
    sout = "K_{"
    indices = []
    for i in range(1, l+1):
        indices.append(int(match_obj.group(i)))
    indices.sort()
    for i in indices:
        sout += str(i)
    return sout + "}"

def combineKNeg(match_obj):
    l = len(match_obj.groups())
    sout = "K_{"
    indices = []
    for i in range(1, l+1):
        indices.append(int(match_obj.group(i)))
    indices.sort()
    for i in indices:
        sout += str(i)
    return sout + "}^{-1}"

def splitE(match_obj):
    return "E_{" + match_obj.group(1) + "} E_{" + match_obj.group(1) + "}"

def splitF(match_obj):
    return "F_{" + match_obj.group(1) + "} F_{" + match_obj.group(1) + "}"

def splitK(match_obj):
    return "K_{" + match_obj.group(1) + "} K_{" + match_obj.group(1) + "}"

def splitKNeg(match_obj):
    return "KNeg_{" + match_obj.group(1) + "} KNeg_{" + match_obj.group(1) + "}"

out = re.sub(r"E\_\{(\d)\}\^\{2\}", splitE, input) # split E
out = re.sub(r"F\_\{(\d)\}\^\{2\}", splitF, out) # split F
out = re.sub(r"K\_\{(\d)\}\^\{2\}", splitK, out) # split K^2
out = re.sub(r"KNeg\_\{(\d)\}\^\{2\}", splitKNeg, out) # split KNeg^2
out = re.sub(r"frac", "tfrac", out) # frac -> tfrac




for currLen in reversed(range(1, 10)):
    regF = "F\_\{(\d)\}" + (currLen - 1) * " F\_\{(\d)\}" 
    regE = "E\_\{(\d)\}" + (currLen - 1) * " E\_\{(\d)\}"
    regK = "K\_\{(\d)\}" + (currLen - 1) * " K\_\{(\d)\}"
    regKNeg = "KNeg\_\{(\d)\}" + (currLen - 1) * " KNeg\_\{(\d)\}"
    out = re.sub(regF, combineF, out) # combine F indices
    out = re.sub(regE, combineE, out) # combine E indices
    out = re.sub(regK, combineK, out) # combine F indices
    out = re.sub(regKNeg, combineKNeg, out) # combine E indices


print(out)



def reverse(match_obj):
    l = len(match_obj.groups())
    str = ""
    for i in range(1, l + 1):
        str += "F"
        str += match_obj.group(l + 1 - i)
        if (i < l):
            str += "*"
    return str

#fDual5112 = re.sub(r"F(\d)\*F(\d)\*F(\d)\*F(\d)\*F(\d)\*F(\d)", reverse, eDual5132)

#print(fDual5112)

