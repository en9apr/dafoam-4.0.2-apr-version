/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v4

\*---------------------------------------------------------------------------*/

#include "DAInputPatchVx.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(DAInputPatchVx, 0);
addToRunTimeSelectionTable(DAInput, DAInputPatchVx, dictionary);
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DAInputPatchVx::DAInputPatchVx(
    const word inputName,
    const word inputType,
    fvMesh& mesh,
    const DAOption& daOption,
    const DAModel& daModel,
    const DAIndex& daIndex)
    : DAInput(
        inputName,
        inputType,
        mesh,
        daOption,
        daModel,
        daIndex)
{
}

void DAInputPatchVx::run(const scalarList& input)
{
    /*
    Description:
        Assign input[0] (Vx) directly to the U field on specified patches.
    */

    // NOTE: we need to first update DAGlobalVar::patchVelocity here, so that daFunction-force
    // can use it to compute the force direction.
    /*DAGlobalVar& globalVar =
        const_cast<DAGlobalVar&>(mesh_.thisDb().lookupObject<DAGlobalVar>("DAGlobalVar"));
    globalVar.patchVelocity[0] = input[0];
    globalVar.patchVelocity[1] = 0.0;
    */
    wordList patchNames;
    dictionary patchVSubDict = daOption_.getAllOptions().subDict("inputInfo").subDict(inputName_);
    patchVSubDict.readEntry<wordList>("patches", patchNames);
/*
#ifndef CODI_ADR
    Info << "DAInputPatchVx. " << endl;
    Info << "DAInputPatchVx: Setting Vx = " << input[0] 
         << " at patches " << patchNames << endl;
#endif
*/
    // the streamwise axis of aoa, aoa = tan( U_normal/U_flow )
       
    
    /*word flowAxis = patchVSubDict.getWord("flowAxis");
    word normalAxis = patchVSubDict.getWord("normalAxis");
    scalar UMag = input[0];

    HashTable<label> axisIndices;
    axisIndices.set("x", 0);
    axisIndices.set("y", 1);
    axisIndices.set("z", 2);
    label flowAxisIndex = axisIndices[flowAxis];
    label normalAxisIndex = axisIndices[normalAxis];
    */
    //label flowAxisIndex = 0;
    
    scalar Vx = input[0];
    
    volVectorField& U = const_cast<volVectorField&>(mesh_.thisDb().lookupObject<volVectorField>("U"));

    //scalar aoaRad = input[1] * constant::mathematical::pi / 180.0;
    //scalar UxNew = UMag * cos(aoaRad);
    
    
    //scalar UyNew = 0.0;
    
    // now we assign UxNew and UyNew to the U patches
    forAll(patchNames, idxI)
    {
        word patchName = patchNames[idxI];
        label patchI = mesh_.boundaryMesh().findPatchID(patchName);
        if (mesh_.boundaryMesh()[patchI].size() > 0)
        {
            if (U.boundaryField()[patchI].type() == "fixedValue")
            {
                forAll(U.boundaryField()[patchI], faceI)
                {
                    U.boundaryFieldRef()[patchI][faceI][0] = Vx;
                    //U.boundaryFieldRef()[patchI][faceI][normalAxisIndex] = UyNew;
                    //U.boundaryFieldRef()[patchI][faceI][flowAxisIndex] = input[0];
                }
            }
            else if (U.boundaryField()[patchI].type() == "inletOutlet")
            {
                mixedFvPatchField<vector>& inletOutletPatch =
                    refCast<mixedFvPatchField<vector>>(U.boundaryFieldRef()[patchI]);

                forAll(U.boundaryField()[patchI], faceI)
                {
                    inletOutletPatch.refValue()[faceI][0] = Vx;
                    //inletOutletPatch.refValue()[faceI][normalAxisIndex] = UyNew;
                    //inletOutletPatch.refValue()[faceI][flowAxisIndex] = input[0];
                }
            }
            else
            {
                FatalErrorIn("DAInputPatchVx::run")
                    << "patch type not valid! only support fixedValue or inletOutlet"
                    << exit(FatalError);
            }
        }
    }
    U.correctBoundaryConditions();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
