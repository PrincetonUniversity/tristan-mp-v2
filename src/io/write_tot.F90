#include "../defs.F90"

module m_writetot
  use m_globalnamespace, only: mpi_rank
  use m_outputnamespace, only: params_enable, params_enable, prtl_tot_enable,&
                             & flds_tot_enable, tot_output_index,&
                             & diag_enable, spectra_enable

  use m_aux, only: printDiag
  use m_outputlogistics, only: prepareOutput
  use m_writetotflds, only: writeFields
  use m_writeparams, only: writeParams
  use m_writetotprtl, only: writeParticles
  use m_writediagnostics, only: writeDiagnostics
  use m_writespectra, only: writeSpectra

  implicit none
contains
  subroutine writeTotOutput(time)
    implicit none
    integer, intent(in)        :: time
    integer                    :: step, ierr

    call prepareOutput()

    step = tot_output_index

    if (params_enable) then
      call writeParams(step, time)
    end if

    if (flds_tot_enable) then
      call writeFields(step, time)
    end if

    if (prtl_tot_enable) then
      call writeParticles(step, time)
    end if

    if (spectra_enable) then
      call writeSpectra(step, time)
    end if

    if (diag_enable) then
      call writeDiagnostics(step, time)
    end if

    call printDiag((mpi_rank .eq. 0), "tot_output()", .true.)
    tot_output_index = tot_output_index + 1
  end subroutine writeTotOutput

end module m_writetot
