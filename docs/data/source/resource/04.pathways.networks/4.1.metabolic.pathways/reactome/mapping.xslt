<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  Reactome -> Unified Metabolic Pathways Schema Mapping
  ============================================================
  Source: ./schema.json (Reactome Content Service API)
  Target: ../schema.json (Unified 4.1 Metabolic Pathways Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────┬──────────────┐
  │ Source Field                │ Target Field            │ Transform    │
  ├─────────────────────────────┼─────────────────────────┼──────────────┤
  │ stId                        │ pathway_id              │ Direct       │
  │ displayName                 │ pathway_name            │ Direct       │
  │ speciesName                 │ organism                │ Direct       │
  │ dbId                        │ reactome_db_id          │ Direct       │
  │ stIdVersion                 │ reactome_st_id_version  │ Direct       │
  │ schemaClass                 │ reactome_schema_class   │ Direct       │
  │ compartment                 │ compartment             │ Array        │
  │ goBiologicalProcess         │ go_biological_process   │ Object       │
  │ literatureReference         │ literature_references   │ Restructure  │
  │ input                       │ input                   │ Array        │
  │ output                      │ output                  │ Array        │
  │ catalystActivity            │ catalyst_activity       │ Array        │
  │ hasEvent                    │ has_event               │ Array        │
  │ referenceEntity             │ reference_entity        │ Object       │
  │ isInDisease                 │ is_in_disease           │ Boolean      │
  │ isInferred                  │ is_inferred             │ Boolean      │
  └─────────────────────────────┴─────────────────────────┴──────────────┘

  NULL HANDLING:
  - Empty strings "" -> null
  - Missing optional fields -> null

  NOTES:
  - Reactome uses stable identifiers (stId) in R-HSA-NNNNN format
  - dbId is the internal database ID (integer)
  - schemaClass indicates entity type (Pathway, Reaction, Complex, etc.)
-->
<xsl:stylesheet version="3.0"
    xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
    xmlns:xs="http://www.w3.org/2001/XMLSchema"
    xmlns:fn="http://www.w3.org/2005/xpath-functions"
    xmlns:map="http://www.w3.org/2005/xpath-functions/map"
    xmlns:array="http://www.w3.org/2005/xpath-functions/array"
    xmlns:local="http://local.functions"
    exclude-result-prefixes="xs fn map array local">

  <xsl:output method="json" indent="yes"/>

  <!-- ============================================================
       Entry Point: Parse JSON input
       ============================================================ -->
  <xsl:template name="xsl:initial-template">
    <xsl:param name="input" as="xs:string"/>
    <xsl:variable name="source" select="json-to-xml($input)"/>
    <xsl:apply-templates select="$source/fn:map"/>
  </xsl:template>

  <!-- ============================================================
       Main Record Transformation
       ============================================================ -->
  <xsl:template match="fn:map">
    <fn:map>
      <!-- ========== DIRECT FIELD MAPPINGS ========== -->

      <!-- Primary identifier (R-HSA-NNNNN) -->
      <fn:string key="pathway_id">
        <xsl:value-of select="fn:string[@key='stId']"/>
      </fn:string>

      <!-- Pathway name -->
      <fn:string key="pathway_name">
        <xsl:value-of select="fn:string[@key='displayName']"/>
      </fn:string>

      <!-- Organism -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='speciesName'] and fn:string[@key='speciesName'] != ''">
          <fn:string key="organism">
            <xsl:value-of select="fn:string[@key='speciesName']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="organism"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== REACTOME-SPECIFIC FIELDS ========== -->

      <!-- Internal database ID -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='dbId']">
          <fn:number key="reactome_db_id">
            <xsl:value-of select="fn:number[@key='dbId']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="reactome_db_id"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Stable ID (same as pathway_id) -->
      <fn:string key="reactome_st_id">
        <xsl:value-of select="fn:string[@key='stId']"/>
      </fn:string>

      <!-- Versioned stable ID -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='stIdVersion'] and fn:string[@key='stIdVersion'] != ''">
          <fn:string key="reactome_st_id_version">
            <xsl:value-of select="fn:string[@key='stIdVersion']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="reactome_st_id_version"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Schema class -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='schemaClass'] and fn:string[@key='schemaClass'] != ''">
          <fn:string key="reactome_schema_class">
            <xsl:value-of select="fn:string[@key='schemaClass']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="reactome_schema_class"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Compartments -->
      <xsl:choose>
        <xsl:when test="fn:array[@key='compartment']">
          <fn:array key="compartment">
            <xsl:for-each select="fn:array[@key='compartment']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="compartment"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- GO Biological Process -->
      <xsl:choose>
        <xsl:when test="fn:map[@key='goBiologicalProcess']">
          <fn:map key="go_biological_process">
            <xsl:copy-of select="fn:map[@key='goBiologicalProcess']/*"/>
          </fn:map>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="go_biological_process"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Literature references -->
      <xsl:choose>
        <xsl:when test="fn:array[@key='literatureReference']">
          <fn:array key="literature_references">
            <xsl:for-each select="fn:array[@key='literatureReference']/fn:map">
              <fn:map>
                <xsl:if test="fn:number[@key='pubMedIdentifier']">
                  <fn:number key="pmid">
                    <xsl:value-of select="fn:number[@key='pubMedIdentifier']"/>
                  </fn:number>
                </xsl:if>
                <xsl:if test="fn:string[@key='title']">
                  <fn:string key="title">
                    <xsl:value-of select="fn:string[@key='title']"/>
                  </fn:string>
                </xsl:if>
              </fn:map>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="literature_references"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Input (substrates) -->
      <xsl:choose>
        <xsl:when test="fn:array[@key='input']">
          <fn:array key="input">
            <xsl:for-each select="fn:array[@key='input']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="input"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Output (products) -->
      <xsl:choose>
        <xsl:when test="fn:array[@key='output']">
          <fn:array key="output">
            <xsl:for-each select="fn:array[@key='output']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="output"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Catalyst activity -->
      <xsl:choose>
        <xsl:when test="fn:array[@key='catalystActivity']">
          <fn:array key="catalyst_activity">
            <xsl:for-each select="fn:array[@key='catalystActivity']/fn:map">
              <fn:map>
                <xsl:copy-of select="./*"/>
              </fn:map>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="catalyst_activity"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Has event (child pathways/reactions) -->
      <xsl:choose>
        <xsl:when test="fn:array[@key='hasEvent']">
          <fn:array key="has_event">
            <xsl:for-each select="fn:array[@key='hasEvent']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="has_event"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Reference entity -->
      <xsl:choose>
        <xsl:when test="fn:map[@key='referenceEntity']">
          <fn:map key="reference_entity">
            <xsl:copy-of select="fn:map[@key='referenceEntity']/*"/>
          </fn:map>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="reference_entity"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Boolean flags -->
      <xsl:choose>
        <xsl:when test="fn:boolean[@key='isInDisease']">
          <fn:boolean key="is_in_disease">
            <xsl:value-of select="fn:boolean[@key='isInDisease']"/>
          </fn:boolean>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="is_in_disease"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:boolean[@key='isInferred']">
          <fn:boolean key="is_inferred">
            <xsl:value-of select="fn:boolean[@key='isInferred']"/>
          </fn:boolean>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="is_inferred"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">Reactome</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:string[@key='stId']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

</xsl:stylesheet>
