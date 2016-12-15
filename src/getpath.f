C
        FUNCTION getpath(f)
              CHARACTER(len=255) :: getpath
              CHARACTER(len=*)   :: f
C
              CHARACTER(len=255) :: DATAPATH
              CALL get_environment_variable("ELSEPA_DATA", DATAPATH)
C
              getpath = TRIM(DATAPATH)//"/"//f
              RETURN
        END FUNCTION
