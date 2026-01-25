<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  JASPAR -> Unified Regulatory Networks Schema Mapping
  ============================================================
  Source: ./schema.json (JASPAR TF binding profiles)
  Target: ../schema.json (Unified 4.6 Regulatory Networks Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────┬──────────────┐
  │ Source Field                │ Target Field            │ Transform    │
  ├─────────────────────────────┼─────────────────────────┼──────────────┤
  │ matrix_id                   │ regulator_id            │ Direct       │
  │ name                        │ regulator_name          │ Direct       │
  │ collection                  │ data_collection         │ Direct       │
  │ tax_group                   │ taxonomic_group         │ Direct       │
  │ class                       │ tf_class                │ Direct       │
  │ family                      │ tf_family               │ Direct       │
  │ consensus                   │ binding_consensus       │ Direct       │
  │ pfm                         │ jaspar_pfm              │ Direct       │
  │ pwm                         │ jaspar_pwm              │ Direct       │
  │ species                     │ species_list            │ Restructure  │
  │ binding_sites               │ jaspar_binding_sites    │ Direct       │
  └─────────────────────────────┴─────────────────────────┴──────────────┘

  NULL HANDLING:
  - Empty strings "" -> null
  - Missing optional fields -> null

  NOTES:
  - Matrix IDs follow MA####.# pattern with version
  - PFM contains raw counts, PWM contains log-odds scores
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

      <!-- Matrix/regulator identifier -->
      <fn:string key="regulator_id">
        <xsl:value-of select="fn:string[@key='matrix_id']"/>
      </fn:string>

      <!-- Transcription factor name -->
      <fn:string key="regulator_name">
        <xsl:value-of select="fn:string[@key='name']"/>
      </fn:string>

      <!-- Collection -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='collection'] and fn:string[@key='collection'] != ''">
          <fn:string key="data_collection">
            <xsl:value-of select="fn:string[@key='collection']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="data_collection"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Taxonomic group -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='tax_group'] and fn:string[@key='tax_group'] != ''">
          <fn:string key="taxonomic_group">
            <xsl:value-of select="fn:string[@key='tax_group']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="taxonomic_group"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- TF structural class -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='class'] and fn:string[@key='class'] != ''">
          <fn:string key="tf_class">
            <xsl:value-of select="fn:string[@key='class']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="tf_class"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- TF family -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='family'] and fn:string[@key='family'] != ''">
          <fn:string key="tf_family">
            <xsl:value-of select="fn:string[@key='family']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="tf_family"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Consensus sequence -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='consensus'] and fn:string[@key='consensus'] != ''">
          <fn:string key="binding_consensus">
            <xsl:value-of select="fn:string[@key='consensus']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="binding_consensus"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== ARRAY TRANSFORMATIONS ========== -->

      <!-- Species list -->
      <xsl:choose>
        <xsl:when test="fn:array[@key='species']">
          <fn:array key="species_list">
            <xsl:for-each select="fn:array[@key='species']/fn:map">
              <fn:map>
                <xsl:if test="fn:string[@key='tax_id']">
                  <fn:string key="taxon_id">
                    <xsl:value-of select="fn:string[@key='tax_id']"/>
                  </fn:string>
                </xsl:if>
                <xsl:if test="fn:string[@key='species']">
                  <fn:string key="name">
                    <xsl:value-of select="fn:string[@key='species']"/>
                  </fn:string>
                </xsl:if>
              </fn:map>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="species_list"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== JASPAR-SPECIFIC FIELDS ========== -->

      <!-- Position Frequency Matrix -->
      <xsl:choose>
        <xsl:when test="fn:map[@key='pfm']">
          <fn:map key="jaspar_pfm">
            <xsl:copy-of select="fn:map[@key='pfm']/*"/>
          </fn:map>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="jaspar_pfm"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Position Weight Matrix -->
      <xsl:choose>
        <xsl:when test="fn:map[@key='pwm']">
          <fn:map key="jaspar_pwm">
            <xsl:copy-of select="fn:map[@key='pwm']/*"/>
          </fn:map>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="jaspar_pwm"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Binding sites -->
      <xsl:choose>
        <xsl:when test="fn:array[@key='binding_sites']">
          <fn:array key="jaspar_binding_sites">
            <xsl:for-each select="fn:array[@key='binding_sites']/fn:map">
              <fn:map>
                <xsl:copy-of select="./*"/>
              </fn:map>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="jaspar_binding_sites"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">JASPAR</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:string[@key='matrix_id']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

  <!-- ============================================================
       Reusable Functions
       ============================================================ -->

  <!-- Convert empty/placeholder to empty string (for null handling) -->
  <xsl:function name="local:nullify" as="xs:string">
    <xsl:param name="value" as="xs:string"/>
    <xsl:choose>
      <xsl:when test="$value = ('-', '.', '', 'NA', 'N/A', 'null')"></xsl:when>
      <xsl:otherwise><xsl:value-of select="$value"/></xsl:otherwise>
    </xsl:choose>
  </xsl:function>

</xsl:stylesheet>
