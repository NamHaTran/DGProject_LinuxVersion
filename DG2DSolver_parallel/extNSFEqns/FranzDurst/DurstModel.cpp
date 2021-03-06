#include "./extNSFEqns/FranzDurst/DurstModel.h"
#include <math.h>
#include "DGMath.h"
#include "VarDeclaration.h"


/* GHI CHU:
 * Cac ham trong source code chinh da modified de include model Durst la:
 * process::NSFEq::calcLocalViscousFlux;
 * std::vector<std::vector<double>>calcGaussViscousTerm;
 *
*/
namespace extNSF_Durst {
    /*Variables*/
    bool enable(false);

    void correctViscousTerms(std::vector<std::vector<double>> &diffTerms, std::vector<double> &U, std::vector<double> &dUx, std::vector<double> &dUy)
    {
        std::vector<std::vector<double>> selfDiffTensor(2, std::vector<double>(4, 0.0));
        selfDiffTensor=extNSF_Durst::calcSelfDiffusionTensor(U, dUx, dUy);

        std::tie(diffTerms[0][0], diffTerms[0][1], diffTerms[0][2], diffTerms[0][3]) = extNSF_Durst::calcSelfDiffusionTerms(selfDiffTensor,U[1]/U[0], U[2]/U[0], 1);
        std::tie(diffTerms[1][0], diffTerms[1][1], diffTerms[1][2], diffTerms[1][3]) = extNSF_Durst::calcSelfDiffusionTerms(selfDiffTensor,U[1]/U[0], U[2]/U[0], 2);
    }

    std::tuple<double, double, double, double> calcSelfDiffusionTerms(std::vector< std::vector<double> > &selfDiffusionTensor, double uVal, double vVal, int dir)
    {
        /*selfDiffusionTensor:
        [mDx    tauXx		tauXy		Qx]
        [mDy    tauYx		tauYy		Qy]*/

        double mDx(selfDiffusionTensor[0][0]),
                mDy(selfDiffusionTensor[1][0]),
                tauXy(selfDiffusionTensor[0][2]),
                tauXx(selfDiffusionTensor[0][1]),
                tauYy(selfDiffusionTensor[1][2]),
                Qx(selfDiffusionTensor[0][3]),
                Qy(selfDiffusionTensor[1][3]);
        double viscTerm1(0.0), viscTerm2(0.0), viscTerm3(0.0), viscTerm4(0.0);

        if (dir==1)
        {
            /*1. Ox direction*/
            viscTerm1 = mDx;
            viscTerm2 = tauXx;
            viscTerm3 = tauXy;
            viscTerm4 = tauXx * uVal + tauXy * vVal + Qx - mDx*0.5*(uVal*uVal + vVal*vVal);
        }
        else if (dir==2)
        {
            /*2. Oy direction*/
            viscTerm1 = mDy;
            viscTerm2 = tauXy;
            viscTerm3 = tauYy;
            viscTerm4 = tauXy * uVal + tauYy * vVal + Qy - mDy*0.5*(uVal*uVal + vVal*vVal);
        }
        return std::make_tuple(viscTerm1, viscTerm2, viscTerm3, viscTerm4);
    }

    std::vector<std::vector<double>> calcSelfDiffusionTensor(std::vector<double> &U, std::vector<double> &dUx, std::vector<double> &dUy)
    {
        /*Output matrix has form:
        [mDx    tauXx		tauXy		Qx]
        [mDy    tauYx		tauYy		Qy]
        */
        std::vector<std::vector<double>> OutputMatrix(2, std::vector<double>(4, 0.0));

        double dux(0.0), duy(0.0), dvx(0.0), dvy(0.0), dEx(0.0), dEy(0.0), rhouVal(0.0), rhovVal, rhoVal(0.0), rhoEVal(0.0), TVal;
        double drhox(0.0), drhoy(0.0),
            drhoux(0.0), drhouy(0.0),
            drhovx(0.0), drhovy(0.0),
            drhoEx(0.0), drhoEy(0.0),
            dTx(0.0), dTy(0.0);
        double mDx(0.0), mDy(0.0);

        double uVal(0.0), vVal(0.0), eVal(0.0);
        int index(0);

        rhoVal = U[0];
        rhouVal = U[1];
        rhovVal = U[2];
        rhoEVal = U[3];

        TVal = math::CalcTFromConsvVar(rhoVal,rhouVal,rhovVal,rhoEVal);
        eVal = material::Cv*TVal;

        uVal = rhouVal / rhoVal;
        vVal = rhovVal / rhoVal;

        drhox = dUx[0];
        drhoy = dUy[0];

        drhoux = dUx[1];
        drhouy = dUy[1];

        drhovx = dUx[2];
        drhovy = dUy[2];

        drhoEx = dUx[3];
        drhoEy = dUy[3];

        dux = math::calcRhouvEDeriv(drhoux, drhox, rhouVal, rhoVal);
        duy = math::calcRhouvEDeriv(drhouy, drhoy, rhouVal, rhoVal);

        dvx = math::calcRhouvEDeriv(drhovx, drhox, rhovVal, rhoVal);
        dvy = math::calcRhouvEDeriv(drhovy, drhoy, rhovVal, rhoVal);

        dEx = math::calcRhouvEDeriv(drhoEx, drhox, rhoEVal, rhoVal);
        dEy = math::calcRhouvEDeriv(drhoEy, drhoy, rhoEVal, rhoVal);

        dTx = math::calcTDeriv(dEx, dux, dvx, uVal, vVal);
        dTy = math::calcTDeriv(dEy, duy, dvy, uVal, vVal);

        /*Self mass diffusion term*/
        //Luu y cac dao ham deu da nhan mu
        mDx = -(drhox/rhoVal + dTx/(2*TVal));
        mDy = -(drhoy/rhoVal + dTy/(2*TVal));
        OutputMatrix[0][0]=mDx;
        OutputMatrix[0][1]=mDy;

        /*calculate stresses*/
        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                index = 10 * (i + 1) + (j + 1);
                if (index == 11)  //x normal stress (tau_xx)
                {
                    OutputMatrix[i][j+1] = extNSF_Durst::calcDiffusionStressComponent(index, mDx*uVal, mDy*vVal);
                }
                else if (index == 22)  //y normal stress (tau_yy)
                {
                    OutputMatrix[i][j+1] = extNSF_Durst::calcDiffusionStressComponent(index, mDy*vVal, mDx*uVal);
                }
                else  //shear stress (tau_xy)
                {
                    OutputMatrix[i][j+1] = extNSF_Durst::calcDiffusionStressComponent(index, mDx*vVal, mDy*uVal);
                }
            }
        }

        /*Sel-diffusion molecular transport of heat*/
        OutputMatrix[0][3] = mDx*eVal;
        OutputMatrix[1][3] = mDy*eVal;
        return OutputMatrix;
    }

    double calcDiffusionStressComponent(int index, double fstTerm, double sndTerm)
    {
        /*Formula of stress
        tau_xx = (4/3)*mDx*ux - (2/3)mDy*uy
        tau_yy = (4/3)*mDy*uy - (2/3)mDx*ux
        tau_xy = mDx*uy + mDy*ux

        input variables:
        index			fstTerm         sndTerm
        11 (tau_xx)		mDx*ux			mDy*Uy
        22 (tau_yy)		mDy*Uy			mDx*ux
        12 (tau_xy)		mDx*uy			mDy*ux

        Dau - trong cac term deu da co trong mDx mDy
        */
        double tau(0.0);
        if (index == 11 || index == 22)
        {
            tau = (4.0 / 3.0)*fstTerm - (2.0 / 3.0)*sndTerm;
        }
        else if (index == 12 || index == 21)
        {
            tau = fstTerm + sndTerm;
        }
        return tau;
    }
}

