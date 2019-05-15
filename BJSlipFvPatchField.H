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

#ifndef BJSlipFvPatchVectorField_H
#define BJSlipFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "transformFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class radiationConvectionFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class BJSlipFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- fluid density
        //scalar rho_;

        //- exponent of non-linear Navier slip condition
        //scalar n_;

        //- slip factor
        scalar slipFactor_;

 	//- relaxation factor
        scalar relaxationFactor_;

	//-  flow direction
	vector n (1 0 0);

	scalar alpha_
	(
		parameter.H.lookup(alpha)
	);

	scalar epsilon
	(
		parameter.H.lookup(epsilon)
	);


public:

    //- Runtime type information
    TypeName("BJSlip");


    // Constructors

        //- Construct from patch and internal field
        BJSlipFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        BJSlipFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given BJSlipFvPatchField
        //  onto a new patch
        BJSlipFvPatchVectorField
        (
            const BJSlipFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

	//- Construct as copy
	BJSlipFvPatchVectorField
	(
	    const BJSlipFvPatchVectorField&
	);

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new BJSlipFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        BJSlipFvPatchVectorField
        (
            const BJSlipFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new BJSlipFvPatchVectorField(*this, iF)
            );
        }


    // Member functions


        //- Return the fluid density
        //const scalar& rho() const
        //{
          //  return rho_;
        //}

        //- Return valuefraction slipFactor
        const scalar& slipFactor() const
        {
            return slipFactor_;
        }

        //- Return n
        //const scalar& n() const
        //{
          //  return n_;
        //}

        //- Return relaxationFactor

        const scalar& relaxationFactor() const
        {
            return relaxationFactor_;
        }

        //- Return non-const acess to the fluid density
       // scalar& rho()
        //{
          //  return rho_;
        //}

        //- Return non-const acess to the valuefraction slipFactor
        scalar& slipFactor()
        {
            return slipFactor_;
        }

        //- Return non-const acess to the n value
       // scalar& n()
        //{
          //  return n_;
        //}

        //- Return non-const acess to the relaxationFactor value

        scalar& relaxationFactor()
        {
            return relaxationFactor_;
        }

        //- Update coefficients
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
