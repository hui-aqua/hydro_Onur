/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
 No.
Class
    Foam::netPanel

Description

SourceFiles
    netPanel.H
    netPanel.C

\*---------------------------------------------------------------------------*/
#ifndef netPanel_H
#define netPanel_H
// class guards that prevents double declaration

#include "dictionary.H"
#include "cellSet.H"
#include "fvMesh.H"
#include "scalar.H"
#include "fvMatrices.H"
#include "volFields.H"
#include "SortableList.H"
// might be used
// #include "pointField.H"
// #include "volFields.H"
// #include "fvMatricesFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    /*---------------------------------------------------------------------------*\
                         Class netPanel Declaration
\*---------------------------------------------------------------------------*/

    class netPanel
    {
    protected:
        //- Private data
        const dictionary netDict_memb;
        const scalar thickness_memb; // thickness of the net panel in CFD, thickness of the porous media.
        const scalar Sn_memb;        // solidity of the netting
        const scalar ML_memb;        // length between neighbor knots
        const scalar dw_memb;        // diameter of the twine in netting
        const scalar fluidrho_memb;

        const scalar updateInterval_memb; //  update interval for velocity to speed the simulation
        const scalar ropeEnhance_memb;    // 0: no enhance, might be a hole;
                                          // 1: the thickness of twine is the same with thickness of net;
                                          // 2: enhanced, the thickness of twine is twice the thickness of net.
        // input
        List<vector> structuralElements_memb;  // screen model of netting
        List<vector> structuralPositions_memb; // positions of each nodes
        List<vector> structuralForces_memb;    // hydrodynamic forces on all the screen elements
        List<scalar> position_range;           // [[x_min,x_max],[y_min,y_max],[z_min,z_max]]
        // output
        List<vector> fluidVelocity_memb; // velocity on the center of all the screen elements

        //-Private member function
        vector calcNorm(
            const point &point_a,
            const point &point_b,
            const point &point_c) const;

        vector calcNorm(
            const point &point_a,
            const point &point_b,
            const point &point_c,
            const vector &fluidVelocity) const;

        scalar calcArea(
            const point &point_a,
            const point &point_b,
            const point &point_c) const;

        scalar calcDist(
            const point &point_a,
            const point &point_b) const;

        bool isInPorousZone(
            const point &x,
            const vector &point_a,
            const vector &point_b,
            const vector &point_c) const;

        bool isInPorous_line(
            const point &x,
            const vector &point_a,
            const vector &point_b,
            const vector &point_c) const;

        scalar calcDistanceFromPoint2Panel(
            const point &x,
            const vector &structuralElementi) const;

    public:
        // Constructors
        netPanel(
            const dictionary &netDict);

        //- Destructor
        ~netPanel();

        // public Member Functions
        void addResistance(
            volScalarField &porosityField,
            fvVectorMatrix &UEqn,
            const fvMesh &mesh) const;

        void updatePoroField(
            volScalarField &porosityField,
            const fvMesh &mesh) const;

        void updateVelocity(
            const List<pointField> &gatheredU,
            const List<pointField> &gatheredmesh,
            const scalar &thresholdLength);

        scalar porosity() const
        {
            return 1 - (0.7853981633974483) * (dw_memb / thickness_memb) * Sn_memb; // 0.785... is pi/4
        };

        scalar Sn() const
        {
            return Sn_memb;
        };

        List<vector> posi() const
        {
            return structuralPositions_memb;
        };

        List<vector> surf() const
        {
            return structuralElements_memb;
        };
        List<vector> FluidU() const
        {
            return fluidVelocity_memb;
        };
        List<vector> Fhout() const
        {
            return structuralForces_memb;
        };

        void readPosi(const dictionary &structuralPositions);

        void readSurf(const dictionary &structuralElements);

        void readForce(
            const scalar &time_foam,
            const dictionary &structuralForces_memb);

        // Friend Functions

        // Friend Operators
    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#include "netPanelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //r&  t
