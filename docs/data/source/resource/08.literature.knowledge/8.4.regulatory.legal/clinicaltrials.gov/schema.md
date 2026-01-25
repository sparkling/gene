---
id: schema-clinicaltrials-gov
title: "ClinicalTrials.gov Schema Documentation"
type: schema
parent: README.md
last_updated: 2026-01-24
status: final
tags: [schema, database, clinical-trials, regulatory]
---

# ClinicalTrials.gov - Schema Documentation

## TL;DR

ClinicalTrials.gov provides comprehensive clinical trial data via a REST API (v2) and bulk downloads. The data model captures study design, eligibility, interventions, outcomes, and results.

## Study Object (API v2)

```json
{
  "protocolSection": {
    "identificationModule": {...},
    "statusModule": {...},
    "sponsorCollaboratorsModule": {...},
    "oversightModule": {...},
    "descriptionModule": {...},
    "conditionsModule": {...},
    "designModule": {...},
    "armsInterventionsModule": {...},
    "outcomesModule": {...},
    "eligibilityModule": {...},
    "contactsLocationsModule": {...},
    "referencesModule": {...}
  },
  "resultsSection": {
    "participantFlowModule": {...},
    "baselineCharacteristicsModule": {...},
    "outcomeMeasuresModule": {...},
    "adverseEventsModule": {...}
  },
  "derivedSection": {
    "miscInfoModule": {...},
    "conditionBrowseModule": {...},
    "interventionBrowseModule": {...}
  },
  "hasResults": true
}
```

## Identification Module

```json
{
  "identificationModule": {
    "nctId": "NCT00000001",
    "orgStudyIdInfo": {
      "id": "SPONSOR-001"
    },
    "secondaryIdInfos": [
      {
        "id": "2020-001234-12",
        "type": "EUDRACT_NUMBER"
      }
    ],
    "organization": {
      "fullName": "Sponsor Name",
      "class": "INDUSTRY"
    },
    "briefTitle": "Short Title of the Study",
    "officialTitle": "Full Official Title of the Study",
    "acronym": "STUDY"
  }
}
```

## Status Module

```json
{
  "statusModule": {
    "statusVerifiedDate": "2024-01",
    "overallStatus": "RECRUITING",
    "expandedAccessInfo": {
      "hasExpandedAccess": false
    },
    "startDateStruct": {
      "date": "2023-06-01",
      "type": "ACTUAL"
    },
    "primaryCompletionDateStruct": {
      "date": "2025-12-31",
      "type": "ESTIMATED"
    },
    "completionDateStruct": {
      "date": "2026-06-30",
      "type": "ESTIMATED"
    },
    "studyFirstSubmitDate": "2023-05-15",
    "studyFirstSubmitQcDate": "2023-05-20",
    "studyFirstPostDateStruct": {
      "date": "2023-05-25",
      "type": "ACTUAL"
    },
    "lastUpdateSubmitDate": "2024-01-10",
    "lastUpdatePostDateStruct": {
      "date": "2024-01-15",
      "type": "ACTUAL"
    }
  }
}
```

### Overall Status Values

| Status | Description |
|--------|-------------|
| NOT_YET_RECRUITING | Approved but not started |
| RECRUITING | Actively enrolling |
| ENROLLING_BY_INVITATION | Invitation only |
| ACTIVE_NOT_RECRUITING | Ongoing, not recruiting |
| SUSPENDED | Temporarily halted |
| TERMINATED | Ended early |
| COMPLETED | Finished |
| WITHDRAWN | Withdrawn before enrollment |
| UNKNOWN | Status not verified |

## Design Module

```json
{
  "designModule": {
    "studyType": "INTERVENTIONAL",
    "phases": ["PHASE3"],
    "designInfo": {
      "allocation": "RANDOMIZED",
      "interventionModel": "PARALLEL",
      "interventionModelDescription": "Two-arm parallel design",
      "primaryPurpose": "TREATMENT",
      "maskingInfo": {
        "masking": "QUADRUPLE",
        "maskingDescription": "Double-blind",
        "whoMasked": ["PARTICIPANT", "CARE_PROVIDER", "INVESTIGATOR", "OUTCOMES_ASSESSOR"]
      }
    },
    "enrollmentInfo": {
      "count": 500,
      "type": "ESTIMATED"
    }
  }
}
```

### Study Types

| Type | Description |
|------|-------------|
| INTERVENTIONAL | Tests intervention |
| OBSERVATIONAL | Observes outcomes |
| EXPANDED_ACCESS | Compassionate use |

### Phases

| Phase | Description |
|-------|-------------|
| EARLY_PHASE1 | Exploratory trials |
| PHASE1 | Safety, dosing |
| PHASE2 | Efficacy signals |
| PHASE3 | Confirmatory |
| PHASE4 | Post-market |
| NA | Not applicable |

## Arms and Interventions

```json
{
  "armsInterventionsModule": {
    "armGroups": [
      {
        "label": "Drug A",
        "type": "EXPERIMENTAL",
        "description": "Drug A 100mg daily",
        "interventionNames": ["Drug: Drug A"]
      },
      {
        "label": "Placebo",
        "type": "PLACEBO_COMPARATOR",
        "description": "Matching placebo daily",
        "interventionNames": ["Drug: Placebo"]
      }
    ],
    "interventions": [
      {
        "type": "DRUG",
        "name": "Drug A",
        "description": "Drug A is a selective inhibitor...",
        "armGroupLabels": ["Drug A"],
        "otherNames": ["Generic Name"]
      }
    ]
  }
}
```

### Intervention Types

| Type | Description |
|------|-------------|
| DRUG | Pharmaceutical |
| DEVICE | Medical device |
| BIOLOGICAL | Biologic/vaccine |
| PROCEDURE | Surgical/medical |
| RADIATION | Radiation therapy |
| BEHAVIORAL | Behavioral intervention |
| GENETIC | Gene therapy |
| DIETARY_SUPPLEMENT | Supplements |
| COMBINATION_PRODUCT | Combined |
| DIAGNOSTIC_TEST | Diagnostic |
| OTHER | Other |

## Eligibility Module

```json
{
  "eligibilityModule": {
    "eligibilityCriteria": "Inclusion Criteria:\\n- Age >= 18\\n- Confirmed diagnosis...\\n\\nExclusion Criteria:\\n- Prior treatment with...\\n- Active infection...",
    "healthyVolunteers": false,
    "sex": "ALL",
    "minimumAge": "18 Years",
    "maximumAge": "75 Years",
    "stdAges": ["ADULT", "OLDER_ADULT"],
    "samplingMethod": "PROBABILITY_SAMPLE",
    "studyPopulation": "Adults with confirmed diagnosis"
  }
}
```

## Outcomes Module

```json
{
  "outcomesModule": {
    "primaryOutcomes": [
      {
        "measure": "Overall Response Rate",
        "description": "Percentage of participants with complete or partial response",
        "timeFrame": "From randomization to 24 weeks"
      }
    ],
    "secondaryOutcomes": [
      {
        "measure": "Progression-Free Survival",
        "description": "Time from randomization to progression or death",
        "timeFrame": "Up to 5 years"
      }
    ]
  }
}
```

## Results Section

```json
{
  "resultsSection": {
    "participantFlowModule": {
      "groups": [
        {
          "id": "FG000",
          "title": "Drug A",
          "description": "Participants receiving Drug A"
        }
      ],
      "periods": [
        {
          "title": "Overall Study",
          "milestones": [
            {
              "type": "STARTED",
              "achievements": [{"groupId": "FG000", "numSubjects": "250"}]
            },
            {
              "type": "COMPLETED",
              "achievements": [{"groupId": "FG000", "numSubjects": "220"}]
            }
          ]
        }
      ]
    },
    "outcomeMeasuresModule": {
      "outcomeMeasures": [
        {
          "title": "Overall Response Rate",
          "type": "PRIMARY",
          "paramType": "NUMBER",
          "unitOfMeasure": "Percentage",
          "classes": [
            {
              "categories": [
                {
                  "measurements": [
                    {"groupId": "OG000", "value": "45.2"}
                  ]
                }
              ]
            }
          ]
        }
      ]
    }
  }
}
```

## Contacts and Locations

```json
{
  "contactsLocationsModule": {
    "centralContacts": [
      {
        "name": "Study Contact",
        "role": "CONTACT",
        "phone": "+1-555-123-4567",
        "email": "study@example.com"
      }
    ],
    "locations": [
      {
        "facility": "Hospital Name",
        "city": "Boston",
        "state": "Massachusetts",
        "zip": "02115",
        "country": "United States",
        "geoPoint": {
          "lat": 42.3601,
          "lon": -71.0589
        },
        "status": "RECRUITING"
      }
    ]
  }
}
```

## See Also

- [Download Documentation](./download.md)
- [ClinicalTrials.gov API](https://clinicaltrials.gov/data-api/api)
