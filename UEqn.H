    surfaceScalarField muEff
    (
        "muEff",
        twoPhaseProperties->muf()
   //   + fvc::interpolate(rho*turbulence->nut())
	      + fvc::interpolate((rho2+pow(alpha1,rebExp)*(rho1-rho2))/(rho2+alpha1*(rho1-rho2))*rho*turbulence->nut())
    );
    // the following is for debugging the correction of muEff according to
    // Reboud
/*    surfaceScalarField muMol
    ( 
        "muMol",
        twoPhaseProperties->muf()
    );
    surfaceScalarField nuTurbf
    ( 
        "nuTurbf",
        fvc::interpolate(turbulence->nut())
    );
    surfaceScalarField rest
    (
        "rest",
	      fvc::interpolate((rho2+pow(alpha1,rebExp)*(rho1-rho2))/(rho2+alpha1*(rho1-rho2))*rho*turbulence->nut())
    );
    surfaceScalarField rhof
    ( 
        "rhof",
        fvc::interpolate(rho)
    );
    surfaceScalarField alpha1f
    ( 
        "alpha1f",
        fvc::interpolate(alpha1)
    );
    forAll(alpha1f.internalField(),cellI)
    {
        if (alpha1f[cellI] < 0.7 && alpha1f[cellI] > 0.3)
        {
            Info << "rebExp = " <<rebExp<<endl; 
            Info << "rho1 = " <<rho1<<endl; 
            Info << "rho2 = " <<rho2<<endl; 
            Info << "rhof = " <<rhof[cellI]<<endl; 
            Info << "alpha1f = " <<alpha1f[cellI]<<endl; 
            Info << "muMol = " <<muMol[cellI]<<endl; 
            Info << "nuTurbf = " <<nuTurbf[cellI]<<endl; 
            Info << "rest = " <<rest[cellI]<<endl; 
            Info << "pow(alpha1f[cellI],rebExp) = " <<pow(alpha1f[cellI],rebExp)<<endl; 
//            Info << "turbulence->nut() = " <<turbulence->nut()<<endl; 
//            Info << "twoPhaseProperties->muf() = " <<twoPhaseProperties->muf()<<endl; 
//            Info << "fvc::interpolate((...*rho*turbulence->nut()) = " <<fvc::interpolate((rho2+pow(alpha1[cellI],rebExp)*(rho1-rho2))/(rho2+alpha1[cellI]*(rho1-rho2))*rho[cellI]*turbulence->nut())<<endl; 
            Info << "cellI = " <<cellI<<endl; 
            Info << "muEff = " <<muEff[cellI]<<endl<<endl; 
        }
    } */

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(rhoPhi, U)
      - fvm::Sp(fvc::ddt(rho) + fvc::div(rhoPhi), U)
      - fvm::laplacian(muEff, U)
      - (fvc::grad(U) & fvc::grad(muEff))
    //- fvc::div(muEff*(fvc::interpolate(dev2(fvc::grad(U))) & mesh.Sf()))
    );

    UEqn.relax();

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1)
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(pd)
                ) * mesh.magSf()
            )
        );
    }
