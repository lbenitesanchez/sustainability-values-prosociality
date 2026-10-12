# Participant data and access

The original participant-level SPSS file is not included. The expected private file is `datos_originales/Base_Final.sav`; its SHA-256 in this audited archive is `8208ee886a79a7348b8cba612c93f275da6e2a393eaf34045e4feb4c78ec9276`. The original 653 questionnaires contain individual responses and demographic information. Authorised users must obtain the unchanged file and the private project dependencies directly from the study authors under the applicable institutional/data-use conditions. This package does not create a new permission to redistribute it. No contact address or public download is invented here.

Raw SPSS CSR1 is Yes=1, No=2. Analytical CSR1_bin is Yes=1, No=0. The analytical dictionary documents institution indicators, sex, study level, age and socioeconomic coding. There is one institution per country; these labels must not be interpreted as nationally representative effects.

Public summary matrices are supporting aggregates, not a substitute for the private input tree needed for ordinal analysis and grouped resampling. Synthetic examples generated with R/00_master.R --mode=synthetic are wholly artificial and do not reproduce published numerical results.
