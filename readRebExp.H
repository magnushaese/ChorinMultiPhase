    IOdictionary RASProperties
    (
        IOobject
        (
            "RASProperties",
            runTime.constant(),
            mesh,
            //IOobject::MUST_READ_IF_MODIFIED,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    
    dimensionedScalar rebExp
    (
		dimensionedScalar::lookupOrDefault
		(
			"reboudExponent",
			RASProperties,
			1.0
		)
	);
	
	Info << "Employing Reboud correction with n = " << rebExp.value() << endl;
