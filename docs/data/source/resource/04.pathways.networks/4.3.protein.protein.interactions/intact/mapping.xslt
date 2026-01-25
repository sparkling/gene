<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  IntAct -> Unified PPI Schema Mapping
  ============================================================
  Source: ./schema.json (IntAct PSI-MI TAB 2.7 format)
  Target: ../schema.json (Unified 4.3 PPI Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────┬──────────────┐
  │ Source Field                │ Target Field            │ Transform    │
  ├─────────────────────────────┼─────────────────────────┼──────────────┤
  │ interactionAc               │ intact_interaction_ac   │ Direct       │
  │ idA                         │ interactor_a_id         │ Direct       │
  │ idB                         │ interactor_b_id         │ Direct       │
  │ altA                        │ gene_symbol_a           │ Extract      │
  │ altB                        │ gene_symbol_b           │ Extract      │
  │ detMethod                   │ detection_method        │ Direct       │
  │ type                        │ interaction_type        │ Direct       │
  │ pubid                       │ publication_ids         │ Split array  │
  │ taxidA                      │ organism_a              │ Extract int  │
  │ taxidB                      │ organism_b              │ Extract int  │
  │ intact-miscore              │ intact_miscore          │ Direct       │
  │ expansion                   │ expansion_method        │ Direct       │
  │ bioRoleA                    │ biological_role_a       │ Direct       │
  │ expRoleA                    │ experimental_role_a     │ Direct       │
  │ expRoleB                    │ experimental_role_b     │ Direct       │
  │ typeA                       │ interactor_type_a       │ Direct       │
  │ xrefA                       │ feature_a               │ Split array  │
  │ stoichiometryA              │ stoichiometry_a         │ Direct       │
  │ negative                    │ negative                │ Boolean      │
  │ hostOrganism                │ host_organism           │ Direct       │
  │ parameters                  │ parameters_interaction  │ Split array  │
  │ checksum                    │ checksum_interaction    │ Direct       │
  └─────────────────────────────┴─────────────────────────┴──────────────┘

  NULL HANDLING:
  - "-" -> null
  - Empty strings "" -> null

  NOTES:
  - IntAct uses PSI-MI controlled vocabulary
  - Identifiers prefixed with database name (e.g., "uniprotkb:P04637")
  - MI scores range from 0 to 1
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
      <!-- ========== REQUIRED FIELDS ========== -->

      <!-- Interactor A ID -->
      <fn:string key="interactor_a_id">
        <xsl:value-of select="fn:string[@key='idA']"/>
      </fn:string>

      <!-- Interactor B ID -->
      <fn:string key="interactor_b_id">
        <xsl:value-of select="fn:string[@key='idB']"/>
      </fn:string>

      <!-- ========== GENE SYMBOLS ========== -->

      <xsl:choose>
        <xsl:when test="fn:string[@key='altA'] and fn:string[@key='altA'] != '' and fn:string[@key='altA'] != '-'">
          <fn:string key="gene_symbol_a">
            <xsl:value-of select="fn:string[@key='altA']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="gene_symbol_a"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:string[@key='altB'] and fn:string[@key='altB'] != '' and fn:string[@key='altB'] != '-'">
          <fn:string key="gene_symbol_b">
            <xsl:value-of select="fn:string[@key='altB']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="gene_symbol_b"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== INTACT-SPECIFIC FIELDS ========== -->

      <!-- Interaction accession -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='interactionAc'] and fn:string[@key='interactionAc'] != ''">
          <fn:string key="intact_interaction_ac">
            <xsl:value-of select="fn:string[@key='interactionAc']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="intact_interaction_ac"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- MI score -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='intact-miscore']">
          <fn:number key="intact_miscore">
            <xsl:value-of select="fn:number[@key='intact-miscore']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="intact_miscore"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== DETECTION METHOD AND TYPE ========== -->

      <xsl:choose>
        <xsl:when test="fn:string[@key='detMethod'] and fn:string[@key='detMethod'] != '' and fn:string[@key='detMethod'] != '-'">
          <fn:string key="detection_method">
            <xsl:value-of select="fn:string[@key='detMethod']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="detection_method"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:string[@key='type'] and fn:string[@key='type'] != '' and fn:string[@key='type'] != '-'">
          <fn:string key="interaction_type">
            <xsl:value-of select="fn:string[@key='type']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="interaction_type"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== EXPANSION AND ROLES ========== -->

      <xsl:choose>
        <xsl:when test="fn:string[@key='expansion'] and fn:string[@key='expansion'] != '' and fn:string[@key='expansion'] != '-'">
          <fn:string key="expansion_method">
            <xsl:value-of select="fn:string[@key='expansion']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="expansion_method"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:string[@key='bioRoleA'] and fn:string[@key='bioRoleA'] != '' and fn:string[@key='bioRoleA'] != '-'">
          <fn:string key="biological_role_a">
            <xsl:value-of select="fn:string[@key='bioRoleA']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="biological_role_a"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:string[@key='expRoleA'] and fn:string[@key='expRoleA'] != '' and fn:string[@key='expRoleA'] != '-'">
          <fn:string key="experimental_role_a">
            <xsl:value-of select="fn:string[@key='expRoleA']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="experimental_role_a"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:string[@key='expRoleB'] and fn:string[@key='expRoleB'] != '' and fn:string[@key='expRoleB'] != '-'">
          <fn:string key="experimental_role_b">
            <xsl:value-of select="fn:string[@key='expRoleB']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="experimental_role_b"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== INTERACTOR TYPE ========== -->

      <xsl:choose>
        <xsl:when test="fn:string[@key='typeA'] and fn:string[@key='typeA'] != '' and fn:string[@key='typeA'] != '-'">
          <fn:string key="interactor_type_a">
            <xsl:value-of select="fn:string[@key='typeA']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="interactor_type_a"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== FEATURES AND STOICHIOMETRY ========== -->

      <!-- Features (pipe-delimited) -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='xrefA'] and fn:string[@key='xrefA'] != '' and fn:string[@key='xrefA'] != '-'">
          <fn:array key="feature_a">
            <xsl:for-each select="tokenize(fn:string[@key='xrefA'], '\|')">
              <fn:string><xsl:value-of select="normalize-space(.)"/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="feature_a"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='stoichiometryA']">
          <fn:number key="stoichiometry_a">
            <xsl:value-of select="fn:number[@key='stoichiometryA']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="stoichiometry_a"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== NEGATIVE AND HOST ========== -->

      <xsl:choose>
        <xsl:when test="fn:boolean[@key='negative']">
          <fn:boolean key="negative">
            <xsl:value-of select="fn:boolean[@key='negative']"/>
          </fn:boolean>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="negative"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:string[@key='hostOrganism'] and fn:string[@key='hostOrganism'] != '' and fn:string[@key='hostOrganism'] != '-'">
          <fn:string key="host_organism">
            <xsl:value-of select="fn:string[@key='hostOrganism']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="host_organism"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== PARAMETERS AND CHECKSUM ========== -->

      <xsl:choose>
        <xsl:when test="fn:string[@key='parameters'] and fn:string[@key='parameters'] != '' and fn:string[@key='parameters'] != '-'">
          <fn:array key="parameters_interaction">
            <xsl:for-each select="tokenize(fn:string[@key='parameters'], '\|')">
              <fn:string><xsl:value-of select="normalize-space(.)"/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="parameters_interaction"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:string[@key='checksum'] and fn:string[@key='checksum'] != '' and fn:string[@key='checksum'] != '-'">
          <fn:string key="checksum_interaction">
            <xsl:value-of select="fn:string[@key='checksum']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="checksum_interaction"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== ORGANISM INFO ========== -->

      <xsl:choose>
        <xsl:when test="fn:number[@key='taxidA']">
          <fn:number key="organism_a">
            <xsl:value-of select="fn:number[@key='taxidA']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="organism_a"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='taxidB']">
          <fn:number key="organism_b">
            <xsl:value-of select="fn:number[@key='taxidB']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="organism_b"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== PUBLICATIONS ========== -->

      <xsl:choose>
        <xsl:when test="fn:string[@key='pubid'] and fn:string[@key='pubid'] != '' and fn:string[@key='pubid'] != '-'">
          <fn:array key="publication_ids">
            <xsl:for-each select="tokenize(fn:string[@key='pubid'], '\|')">
              <fn:string><xsl:value-of select="normalize-space(.)"/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="publication_ids"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">IntAct</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:string[@key='interactionAc']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

</xsl:stylesheet>
