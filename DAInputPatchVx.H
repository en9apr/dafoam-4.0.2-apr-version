/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v4

    Description:
        Child class for velocity magnitude and angle of attach at a given patch

\*---------------------------------------------------------------------------*/

#ifndef DAInputPatchVx_H
#define DAInputPatchVx_H

#include "DAInput.H"
#include "addToRunTimeSelectionTable.H"
#include "mixedFvPatchFields.H"
#include "DAGlobalVar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class DAInputPatchVx Declaration
\*---------------------------------------------------------------------------*/

class DAInputPatchVx
    : public DAInput
{

protected:
public:
    TypeName("patchVx");
    // Constructors

    //- Construct from components
    DAInputPatchVx(
        const word inputName,
        const word inputType,
        fvMesh& mesh,
        const DAOption& daOption,
        const DAModel& daModel,
        const DAIndex& daIndex);

    //- Destructor
    virtual ~DAInputPatchVx()
    {
    }

    virtual void run(const scalarList& input);

    virtual label size()
    {
        // [UMag, AoA]
        //return 2;
        return 1;
    }

    virtual label distributed()
    {
        return 0; // not distrubted across processes
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
