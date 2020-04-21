//NOTE: 更新通量
{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    surfaceScalarField phic(mag(phi/mesh.magSf()));
    //surfaceScalarField phic(mixture.cAlpha()*mag(phi/mesh.magSf()));

    // Add the optional isotropic compression contribution
    // if (icAlpha > 0)
    // {
    //     phic *= (1.0 - icAlpha);
    //     phic += (mixture.cAlpha()*icAlpha)*fvc::interpolate(mag(U));
    // }

    // Do not compress interface at non-coupled boundary faces
    // (inlets, outlets etc.)
    surfaceScalarField::Boundary& phicBf = phic.boundaryFieldRef();
    
    forAll(phic.boundaryField(), patchi)
    {
        fvsPatchScalarField& phicp = phicBf[patchi];

        if (!phicp.coupled())
        {
            phicp == 0;
        }
    }

    tmp<surfaceScalarField> phiCN(phi);
//2.3
    //  phic = min(mixture.cAlpha()*phic, max(phic));
//  surfaceScalarField phir(phic*interface.nHatf());
    surfaceScalarField phir(phic*nVecf);

    surfaceScalarField phiH
    (
        fvc::flux
        (
            phiCN(),
            H,//alpha1,
            alphaScheme
        )
      + fvc::flux
        (
            -fvc::flux(-phir, scalar(1) - H, alpharScheme),
            H,//alpha1, 
            alpharScheme
        )
    );


    //MULES::explicitSolve(alpha1, phi, phiAlpha, 1, 0);
   /*  MULES::explicitSolve(rho, H, phi, phiH, 1.0, 0.0); */
    #include "alphaSuSp.H"
    MULES::explicitSolve
    (
        geometricOneField(),
        H,
        phiCN,
        phiH,
        Sp,
        (Su + divU*min(alpha1(), scalar(1)))(),
        oneField(),
        zeroField()
    );
    
    rhoPhiH = phiH*(rho1 - rho2) + phi*rho2;
}
