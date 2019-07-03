# validation

### A wish-list of infromation required for ASKAP spectral line data validation

| Basic Info              | Status | Required file                    | Code variable  | Suggested by |
|-------------------------|:------:|----------------------------------|----------------|--------------|
| SBID                    | ok     | None, input parameter            | sbid           | BQF          |
| Obs-date start          | ok     | metadata/mslist\*txt              | start_obs_date | BQF          |
| Obs-date end            | ok     | metadata/mslist\*txt              | end_obs_date   | BQF          |
| Obs duration (hrs)      | ok     | metadata/mslist\*txt              | duration_hrs   | BQF          |
| ASKAPsoft version       |        | Unknown (fits header history)     | askapsoft      | BQF          |
| Field                   | ok     | metadata/mslist\*txt              | field          | BQF          |
| R.A                     | ok     | metadata/mslist\*txt              | ra             | BQF          |
| Dec                     | ok     | metadata/mslist\*txt              | dec            | BQF          |
| Observed bandwidth      | ok     | metadata/mslist\*txt              | bw             | BQF          |
| Central frequency       | ok     | metadata/mslist\*txt              | cfreq          | BQF          |
| Number of antenna       | ok     | metadata/mslist\*txt              | number_ant     | BQF          |
| Processed channel range | ok     | slurmOutput/\<latest\_executed\>.sh | chan_range     | BQF          |
| Synthesized beam (bmaj) |        | Unknown (fits header)             | bmaj             | BQF          |
| Synthesized beam (bmin) |        | Unknown (fits header)             | bmin               | BQF          
### Validation metrics required for ASKAP pectral line data

| Validation metrics      | Status | Required file                    | Edited by |
|-------------------------|:------:|----------------------------------|--------------|
| Flagging fringe rotation|        | Measurement set                  | BQF          |
| Flagging RFI            |        | Measurement set                  | BQF          |  
| Flagging bandwidth edges|        | Measurement set                  | BQF          |
| Flagging autocorrelations |        | Measurement set                | BQF          |
| Max flux density        |        | CubeStats\*txt                   | BQF          |
| Continuum subtraction   |        |                                  | BQF          |
| Cleaning                |        |                                  | BQF          |
| Sidelobes               |        |                                  | BQF          |
| 1 percentie noise level |        |                                  | BQF          |

Further information: https://jira.csiro.au/browse/ACES-368
Useful resources from SPACS: http://spacs.pbworks.com/w/page/126067640/dataquality
