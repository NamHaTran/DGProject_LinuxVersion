#include "limiterController.h"
#include "VarDeclaration.h"
#include "dynamicVarDeclaration.h"
#include "positivityPreserving.h"
#include "massDiffusion.h"
#include "pAdaptive.h"

namespace limiter
{
    //ham limiter chay sau khi hoan thanh 1 step TVDRK3
    void limiter_1InnerStep()
    {
        if (mathVar::orderElem != 0)
        {
            limiter::pAdaptive::limiter();

            limiter::positivityPreserving::limiter();
        }
    }

    void limiter_1OutterStep()
    {
        if (limitVal::massDiffusion)
        {
            limiter::massDiffusion::limiter();
        }
    }

    /**
     * @brief Function limits Rho by applying Positivity Preserving Limiter.
     */
    void limitRho_PositivityPreserving()
    {
        if (mathVar::orderElem!=0)
        {
            if (limitVal::PositivityPreserving)
            {
                for (int nelem = 0; nelem < meshVar::nelem2D; nelem++)
                {
                    theta1Arr[nelem] = limiter::positivityPreserving::calcTheta1(nelem);
                }
            }
        }
    }
}
