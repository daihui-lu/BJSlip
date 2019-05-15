/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "BJSlipFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "symmTransformField.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"
#include "fieldTypes.H"
#include "incompressible/turbulenceModel/turbulenceModel.H"
#include "RASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

BJSlipFvPatchVectorField::BJSlipFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    //rho_(0),
   // n_(0),
    slipFactor_(0),
    relaxationFactor_(0)
{}


BJSlipFvPatchVectorField::BJSlipFvPatchVectorField
(
    const BJSlipFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
   // rho_(ptf.rho_),
    //n_(ptf.n_),
    slipFactor_(ptf.slipFactor_),
    relaxationFactor_(ptf.relaxationFactor_)
{}


BJSlipFvPatchVectorField::BJSlipFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF),
   // rho_(readScalar(dict.lookup("rho"))),
    //n_(readScalar(dict.lookup("n"))),
    slipFactor_(readScalar(dict.lookup("slipFactor"))),
    relaxationFactor_(readScalar(dict.lookup("relaxationFactor")))
{
    if (dict.found("value"))
    {
        fvPatchField<vector>::operator=
        (
            vectorField("value", dict, p.size())
        );
    }
    else
    {
        // Evaluate the wall velocity
        updateCoeffs();
    }
}


BJSlipFvPatchVectorField::BJSlipFvPatchVectorField
(
    const BJSlipFvPatchVectorField& fcvpvf
)
:
    fixedValueFvPatchField<vector>(fcvpvf),
    //rho_(fcvpvf.rho_),
    //n_(fcvpvf.n_),
    slipFactor_(fcvpvf.slipFactor_),
    relaxationFactor_(fcvpvf.relaxationFactor_)
{}


BJSlipFvPatchVectorField::BJSlipFvPatchVectorField
(
    const BJSlipFvPatchVectorField& fcvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(fcvpvf, iF),
    //rho_(fcvpvf.rho_),
    //n_(fcvpvf.n_),
    slipFactor_(fcvpvf.slipFactor_),
    relaxationFactor_(fcvpvf.relaxationFactor_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void BJSlipFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    //face normal vector
    vectorField nHat = this->patch().nf();  

    // Du/dn
    vectorField  gradient = this->snGrad();
    
    // only tangential components
    gradient = transform(I - sqr(nHat), gradient);
    
    // gradient Direction (since the pow function doesn't take vectors, we later have to multiply the magnitude of the gradient with the gradient direction)
    vectorField gradientDirection = gradient / (mag(gradient) + SMALL);

    // slip velocity of the last iteration
    vectorField u_wallslip_lastIteration = (*this);
    vectorField v_wallslip_lastIteration = (*this);

    const label patchI = patch().index();
    scalarField nuw = 1e-6*mag(patch().nf());  
    
    //if (db().found("turbulenceModel"))
    //{
	//const incompressible::turbulenceModel& turbModel =
          //  db().lookupObject<incompressible::turbulenceModel>
            //(
              //  "turbulenceModel"
            //);

	//Info<< "\nTurbulence found\n" << endl;
	//nuw = rho_*turbModel.nu()().boundaryField()[patchI];
    //}
    //else
    //{
      //  nuw = mag(patch().nf());  
	//Info<< "\nTurbulence NOT found\n" << endl;
    //}
    //Info<<"\n nuw_val = "<< nuw <<endl;    

    // slip velocity in the current iteration
    vectorField u_wallslip = slipFactor_*mag(gradient_face)*gradientDirection_face;
    vectorField v_wallslip = 4*Foam::pow(1+alpha_/epsilon*cellC&n, 2) * alpha_/epsilon/(6*slipFactor_ + 2*Foam::pow(1+alpha_/epsilon*cellC&n, 2));

    //Info <<"\n u_wallslip = "<< u_wallslip << endl;

    //Calculate and set u_wallslip       
     vectorField::operator=(relaxationFactor_*u_wallslip_lastIteration+(1.0-relaxationFactor_)*u_wallslip);
     Info << u_wallslip[50] << endl;
     fixedValueFvPatchVectorField::updateCoeffs();    
}


// Write
void BJSlipFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
   // os.writeKeyword("rho") << rho_ << token::END_STATEMENT << nl;
   // os.writeKeyword("n") << n_ << token::END_STATEMENT << nl;
    os.writeKeyword("slipFactor") << slipFactor_ << token::END_STATEMENT << nl;
    os.writeKeyword("relaxationFactor") << relaxationFactor_ << token::END_STATEMENT << nl;  
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, BJSlipFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
