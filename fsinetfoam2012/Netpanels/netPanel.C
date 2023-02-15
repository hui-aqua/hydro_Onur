/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
 
 No.

\*---------------------------------------------------------------------------*/

#include "netPanel.H"
// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::netPanel::calcArea(
    const point &point_a,
    const point &point_b,
    const point &point_c) const
{
    const vector ab(point_a - point_b);
    const vector ac(point_a - point_c);
    return 0.5 * mag(ac ^ ab);
}

Foam::scalar Foam::netPanel::calcDistanceFromPoint2Panel(
    const point &x,
    const vector &structuralElementi) const
{
    const point point_a = structuralPositions_memb[structuralElementi[0]];
    const point point_b = structuralPositions_memb[structuralElementi[1]];
    const point point_c = structuralPositions_memb[structuralElementi[2]];
    vector panelNorm = calcNorm(point_a, point_b, point_c); // a unit vector to indicate the normal
    scalar dis(mag((x - point_a) & panelNorm));
    //    Info << "The distance from point to net panel is "<< dis << " m." << endl;
    return dis;
}

bool Foam::netPanel::isInPorous_line(
    const point &x,
    const vector &point_a,
    const vector &point_b,
    const vector &point_c) const
{
    bool result(false);

    // Method to enhance the mesh line ()
    const Foam::scalar line1(mag(point_a - point_b));
    const Foam::scalar line2(mag(point_a - point_c));
    const Foam::scalar line3(mag(point_b - point_c));
    Foam::scalar d1(0);
    Foam::scalar d2(0);
    Foam::scalar dis2line(0);
    if (line1 > line2 and line1 > line3)
    {
        // based on line 2
        d1 = (mag(point_a - x));
        d2 = (mag(point_c - x));
        dis2line = mag(calcArea(point_a, point_c, x) * 2.0 / (line2 + Foam::SMALL));
        if (dis2line <= thickness_memb * 0.5 * ropeEnhance_memb and max(d1, d2) <= (thickness_memb * 0.5 * ropeEnhance_memb + line2))
        {
            result = true;
        }
        else
        {
            // based on line3
            d1 = (mag(point_b - x));
            d2 = (mag(point_c - x));
            dis2line = (mag(calcArea(point_b, point_c, x) * 2.0 / (line3 + Foam::SMALL)));
            if (dis2line <= thickness_memb * 0.5 * ropeEnhance_memb and max(d1, d2) <= (thickness_memb * 0.5 * ropeEnhance_memb + line3))
            {
                result = true;
            }
        }
    }
    else
    {
        if (line2 < line3)
        {
            //line 2
            d1 = (mag(point_a - x));
            d2 = (mag(point_c - x));
            dis2line = (mag(calcArea(point_a, point_c, x) * 2.0 / (line2 + Foam::SMALL)));
            if (dis2line <= thickness_memb * 0.5 * ropeEnhance_memb and max(d1, d2) <= (thickness_memb * 0.5 * ropeEnhance_memb + line2))
            {
                result = true;
            }
            else
            {
                //line 1
                d1 = (mag(point_a - x));
                d2 = (mag(point_b - x));
                dis2line = (mag(calcArea(point_a, point_b, x) * 2.0 / (line1 + Foam::SMALL)));
                if (dis2line <= thickness_memb * 0.5 * ropeEnhance_memb and max(d1, d2) <= (thickness_memb * 0.5 * ropeEnhance_memb + line1))
                {
                    result = true;
                }
            }
        }
        else
        {
            //line 3
            d1 = (mag(point_b - x));
            d2 = (mag(point_c - x));
            dis2line = (mag(calcArea(point_b, point_c, x) * 2.0 / (line3 + Foam::SMALL)));
            if (dis2line <= thickness_memb * 0.5 * ropeEnhance_memb and max(d1, d2) <= (thickness_memb * 0.5 * ropeEnhance_memb + line3))
            {
                result = true;
            }
            else
            {
                //line 1
                d1 = (mag(point_a - x));
                d2 = (mag(point_b - x));
                dis2line = (mag(calcArea(point_a, point_b, x) * 2.0 / (line1 + Foam::SMALL)));
                if (dis2line <= thickness_memb * 0.5 * ropeEnhance_memb and max(d1, d2) <= (thickness_memb * 0.5 * ropeEnhance_memb + line1))
                {
                    result = true;
                }
            }
        }
    }
    return result;
}

bool Foam::netPanel::isInPorousZone(
    const point &x,
    const vector &point_a,
    const vector &point_b,
    const vector &point_c) const
{
    bool result(false); // initial value

    vector panelNorm = calcNorm(point_a, point_b, point_c); // a unit vector to indicate the normal
    scalar dis(mag((x - point_a) & panelNorm));
    // define a const scalar as the distance between point x to net panel
    if (dis <= thickness_memb * 0.5) // distance is less than half thickness
    {
        scalar panelarea(calcArea(point_a, point_b, point_c));
        vector projectedPoint(0, 0, 0);        // initial the projected point is 0,0,0
        if (((x - point_a) & panelNorm) < 0.0) // on the side of normal vector
        {
            projectedPoint = (x + panelNorm * dis);
        }
        else
        {
            projectedPoint = (x - panelNorm * dis);
        }
        // projectedPiont is the projected point on the net panel
        scalar panelarea3(calcArea(point_a, point_b, projectedPoint) +
                          calcArea(point_a, projectedPoint, point_c) +
                          calcArea(projectedPoint, point_b, point_c)); //  the area of the three trigular shapes.
        if (panelarea3 <= Foam::SMALL + panelarea)
        {
            result = true;
        }
    }

    return result;
}
Foam::vector Foam::netPanel::calcNorm(
    const point &point_a,
    const point &point_b,
    const point &point_c,
    const vector &fluidVelocity) const
{
    const vector a(point_a - point_b);
    const vector b(point_a - point_c);
    vector norm(a ^ b / (mag(a ^ b) + SMALL));
    if ((norm & fluidVelocity) < 0) // assume the velocity is x+
    {
        norm = -norm;
    }
    return norm; // must be in the same direction with fulow
}

Foam::vector Foam::netPanel::calcNorm(
    const point &point_a,
    const point &point_b,
    const point &point_c) const
{
    const vector a(point_a - point_b);
    const vector b(point_a - point_c);
    const vector norm(a ^ b / (mag(a ^ b) + SMALL));
    return norm; // can point out or in the panel.
}

Foam::scalar Foam::netPanel::calcDist(
    const point &point_a,
    const point &point_b) const
{
    return mag(point_a - point_b);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::netPanel::netPanel(
    const dictionary &netDict)
    : // initial components
      netDict_memb(netDict),
      Sn_memb(readScalar(netDict_memb.subDict("NetInfo1").lookup("Sn"))),
      thickness_memb(readScalar(netDict_memb.subDict("NetInfo1").lookup("PorousMediaThickness"))),
      ML_memb(readScalar(netDict_memb.subDict("NetInfo1").lookup("halfMeshSize"))),
      fluidrho_memb(readScalar(netDict_memb.subDict("NetInfo1").lookup("fluidDensity"))),
      dw_memb(readScalar(netDict_memb.subDict("NetInfo1").lookup("twineDiameter"))),
      updateInterval_memb(readScalar(netDict_memb.subDict("NetInfo1").lookup("velocityUpdateInterval"))),
      ropeEnhance_memb(readScalar(netDict_memb.subDict("NetInfo1").lookup("ropeEnhance")))
{
    // creat the netpanel object
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
Foam::netPanel::~netPanel()
{
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::netPanel::addResistance(
    volScalarField &porosityField,
    fvVectorMatrix &UEqn,
    const fvMesh &mesh) const
{
    Info << ">>> addResistance" << endl;
    forAll(mesh.C(), cellI)
    {
        porosityField[cellI] = 1.0;
    }
    // get the center of all the cells
    const vectorField &centres(mesh.C());
    const scalarField V = mesh.V();
    vectorField &Usource = UEqn.source();
    //  Info << "In addResistance, number of mesh is " << centres.size() << endl;
    // Info << "The structural elements are " << structuralElements_memb << endl;
    // vector sourceforce = structuralForces_memb;
    forAll(centres, cellI)
    {
        forAll(structuralForces_memb, Elementi)
        {
            point p0(structuralPositions_memb[structuralElements_memb[Elementi][0]]);
            point p1(structuralPositions_memb[structuralElements_memb[Elementi][1]]);
            point p2(structuralPositions_memb[structuralElements_memb[Elementi][2]]);
            scalar area(calcArea(p0, p1, p2));

            if (
                isInPorousZone(centres[cellI], p0, p1, p2) || isInPorous_line(centres[cellI], p0, p1, p2))
            {
                Usource[cellI] -= structuralForces_memb[Elementi] * V[cellI];
                porosityField[cellI] = Sn_memb;
                // Info << ">>> addResistance   >>> add source is "<<Usource[cellI] << endl;
            }
        }
    }
}

void Foam::netPanel::updatePoroField(
    volScalarField &porosityField,
    const fvMesh &mesh) const
{
    //    Info << "In updatePoroField, number of mesh is " << (mesh.C()).size() << endl;
    // Info << "The structural elements are " << structuralElements_memb << endl;
    // step1 set all the cell as 1
    Info << ">>> UpdatePorousField" << endl;
    forAll(mesh.C(), cellI)
    {
        porosityField[cellI] = 1.0;
    }
    // get the center of all the cells
    const vectorField &centres(mesh.C());
    //- step 2 assign sn to the proper mesh
    forAll(structuralElements_memb, Elementi) // loop through all the structural emlements
    {
        point p0(structuralPositions_memb[structuralElements_memb[Elementi][0]]);
        point p1(structuralPositions_memb[structuralElements_memb[Elementi][1]]);
        point p2(structuralPositions_memb[structuralElements_memb[Elementi][2]]);

        forAll(centres, cellI) // loop through all the cell,
        {
            if (
                isInPorousZone(centres[cellI], p0, p1, p2) || isInPorous_line(centres[cellI], p0, p1, p2))
            {
                porosityField[cellI] = Sn_memb;
            }
        }
    }
}

// ge the velocity at the net panel center
void Foam::netPanel::updateVelocity(
    const List<pointField> &gatheredU,
    const List<pointField> &gathered_mesh,
    const scalar &thresholdLength)
{
    Info << ">>> UpdateVelocity" << endl;
    List<vector> fluidVelocities(structuralPositions_memb.size(), vector::zero);
    //    const vectorField &centres(mesh.C());
    //    const label &nCells = mesh.nCells();
    //    Info << "In updateVelocity, number of mesh is " << nCells << endl;
    // Info << "The structural elements are " << structuralElements_memb << endl;
    scalar maxDistance(10);                 //started from 2 m ML_memb
    forAll(structuralPositions_memb, Elemi) // loop through all the structural emlements
    {
        maxDistance = thresholdLength * 2; //started from 2 m ML_memb
        vector nearestCell(0.25, 0, 0);
        scalar loops(0);
        forAll(gathered_mesh, processorI) // loop through all the cell,
        {
            if (maxDistance < thresholdLength)
            {
                break;
            }
            forAll(gathered_mesh[processorI], PointI)
            {
                scalar k1(calcDist(gathered_mesh[processorI][PointI], structuralPositions_memb[Elemi]));
                if (k1 < maxDistance)
                {
                    maxDistance = k1;
                    fluidVelocities[Elemi] = gatheredU[processorI][PointI];
                    nearestCell = gathered_mesh[processorI][PointI];
                    loops += 1;
                }
                if (maxDistance < thresholdLength)
                {
                    break;
                }
            }
        }
        //        Info << "After " << loops << " times of loop, the nearest cell is " << nearestCell << "to point " << EP_center <<", and the velocity is "<<fluidVelocities[Elemi]<< "\n"
        //             << endl;

    }
    fluidVelocity_memb = fluidVelocities; // only assige onece
    // Info << "the velocity on elements are  " << fluidVelocity_memb << endl;
}

// * * * * * * * * * * * * * * Communication Functions  * * * * * * * * * * * * * * //
// - the following function is used to communicate with FE solver.
// - initial the structural element
void Foam::netPanel::readSurf(
    const dictionary &structuralElements)
{
    Info << ">>> readSurf" << endl;
    scalar listLength(readScalar(structuralElements.lookup("numOfSurf")));
    List<vector> surf(listLength, vector::zero);
    forAll(surf, i)
    {
        word surf_name("e" + Foam::name(i));
        surf[i] = vector(structuralElements.lookup(surf_name));
    }
    structuralElements_memb = surf;
    // Info << "The structural elements are " << structuralElements_memb << endl;
}

//- update the position and forces
void Foam::netPanel::readPosi(
    const dictionary &structuralPositions)
{
    Info << ">>> readPosi" << endl;
    scalar listLength(readScalar(structuralPositions.lookup("numOfPoint")));
    List<vector> posi(listLength, vector::zero);
    List<scalar> x_array(listLength, 0);
    List<scalar> y_array(listLength, 0);
    List<scalar> z_array(listLength, 0);
    forAll(posi, i)
    {
        word point_name("p" + Foam::name(i));
        posi[i] = vector(structuralPositions.lookup(point_name));
        x_array[i] = posi[i][0];
        y_array[i] = posi[i][1];
        z_array[i] = posi[i][2];
    }
    structuralPositions_memb = posi;

    Info << "X min is " << min(x_array - thickness_memb) << endl;
    Info << "X max is " << max(x_array + thickness_memb) << endl;

    Info << "Y min is " << min(y_array - thickness_memb) << endl;
    Info << "Y max is " << max(y_array + thickness_memb) << endl;

    Info << "Z min is " << min(z_array - thickness_memb) << endl;
    Info << "Z max is " << max(z_array + thickness_memb) << endl;
    tensor position_range(min(x_array - thickness_memb),0,max(x_array + thickness_memb),
                          min(y_array - thickness_memb),0,max(y_array + thickness_memb),
                          min(z_array - thickness_memb),0,max(z_array + thickness_memb));
    Info << "position range  is " << position_range << endl;
}

void Foam::netPanel::readForce(
    const scalar &time_foam,
    const dictionary &structuralForces)
{
    Info << ">>> readForce" << endl;
    scalar listLength(readScalar(structuralForces.lookup("numOfFh")));
    scalar time_FE(readScalar(structuralForces.lookup("timeInFE")));
    if (mag(time_FE - time_foam) > 10)
    {
        Info << "Warning!!! The difference of time in FE and FV solvers exceeds 10 s!\n"
             << endl;
    }

    List<vector> Fh(listLength, vector::zero);
    forAll(Fh, i)
    {
        word force_name("fh" + Foam::name(i));
        Fh[i] = vector(structuralForces.lookup(force_name));
    }
    structuralForces_memb = Fh;
}

// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

// ************************************************************************* //
