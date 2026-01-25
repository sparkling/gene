<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  Dr. Duke's Phytochemical Database -> Natural Products Unified Schema
  ============================================================
  Source: ./schema.md (Dr. Duke's Phytochemical and Ethnobotanical Databases)
  Target: ../schema.json (2.1 Natural Products Unified Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────┬──────────────┐
  │ Source Field                │ Target Field            │ Transform    │
  ├─────────────────────────────┼─────────────────────────┼──────────────┤
  │ common_name                 │ name                    │ Direct       │
  │ cas_number                  │ cas_number              │ Direct       │
  │ biological_activity         │ biological_activities   │ Split array  │
  │ ethnobotanical_use          │ ethnobotanical_uses     │ Object array │
  │ plant_part                  │ plant_part              │ Direct       │
  │ low/high_concentration      │ concentration_range     │ Combine obj  │
  │ plant_name                  │ organism_sources        │ Object array │
  │ molecular_formula           │ molecular_formula       │ Direct       │
  │ molecular_weight            │ molecular_weight        │ Direct       │
  └─────────────────────────────┴─────────────────────────┴──────────────┘

  NULL HANDLING:
  - common_name: "" → null
  - cas_number: "---" or "" → null
  - concentration values: -1 → null

  NOTES:
  - Dr. Duke's focuses on ethnobotanical uses and biological activities
  - Biological activities are extensive (1,900+ types documented)
  - Plant parts specify where in plant compound is found
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
      <!-- ========== BASIC PROPERTIES ========== -->

      <!-- Compound name -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='common_name'])">
          <fn:null key="name"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="name">
            <xsl:value-of select="fn:string[@key='common_name']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- CAS Number -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='cas_number']) or fn:string[@key='cas_number'] = '---'">
          <fn:null key="cas_number"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="cas_number">
            <xsl:value-of select="fn:string[@key='cas_number']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== MOLECULAR PROPERTIES ========== -->

      <!-- Molecular formula -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='molecular_formula'])">
          <fn:null key="molecular_formula"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="molecular_formula">
            <xsl:value-of select="fn:string[@key='molecular_formula']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Molecular weight -->
      <xsl:choose>
        <xsl:when test="not(fn:number[@key='molecular_weight']) or fn:number[@key='molecular_weight'] &lt;= 0">
          <fn:null key="molecular_weight"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:number key="molecular_weight">
            <xsl:value-of select="fn:number[@key='molecular_weight']"/>
          </fn:number>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== BIOLOGICAL ACTIVITIES ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='biological_activity']/fn:string">
          <fn:array key="biological_activities">
            <xsl:for-each select="fn:array[@key='biological_activity']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:when test="fn:string[@key='biological_activity'] and not(local:is-null(fn:string[@key='biological_activity']))">
          <!-- Handle comma-separated string -->
          <fn:array key="biological_activities">
            <xsl:for-each select="tokenize(fn:string[@key='biological_activity'], ';')">
              <fn:string><xsl:value-of select="normalize-space(.)"/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="biological_activities"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== ETHNOBOTANICAL USES ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='ethnobotanical_use']/fn:map">
          <fn:array key="ethnobotanical_uses">
            <xsl:for-each select="fn:array[@key='ethnobotanical_use']/fn:map">
              <fn:map>
                <fn:string key="use">
                  <xsl:value-of select="fn:string[@key='use']"/>
                </fn:string>
                <xsl:if test="fn:string[@key='culture']">
                  <fn:string key="culture">
                    <xsl:value-of select="fn:string[@key='culture']"/>
                  </fn:string>
                </xsl:if>
                <xsl:if test="fn:string[@key='region']">
                  <fn:string key="region">
                    <xsl:value-of select="fn:string[@key='region']"/>
                  </fn:string>
                </xsl:if>
              </fn:map>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="ethnobotanical_uses"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== PLANT PART ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='plant_part'])">
          <fn:null key="plant_part"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="plant_part">
            <xsl:value-of select="fn:string[@key='plant_part']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== CONCENTRATION RANGE ========== -->

      <xsl:choose>
        <xsl:when test="fn:number[@key='low_concentration'] or fn:number[@key='high_concentration']">
          <fn:map key="concentration_range">
            <xsl:if test="fn:number[@key='low_concentration'] and fn:number[@key='low_concentration'] >= 0">
              <fn:number key="low">
                <xsl:value-of select="fn:number[@key='low_concentration']"/>
              </fn:number>
            </xsl:if>
            <xsl:if test="fn:number[@key='high_concentration'] and fn:number[@key='high_concentration'] >= 0">
              <fn:number key="high">
                <xsl:value-of select="fn:number[@key='high_concentration']"/>
              </fn:number>
            </xsl:if>
            <fn:string key="unit">ppm</fn:string>
          </fn:map>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="concentration_range"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== ORGANISM SOURCES ========== -->

      <xsl:choose>
        <xsl:when test="fn:string[@key='plant_name'] and not(local:is-null(fn:string[@key='plant_name']))">
          <fn:array key="organism_sources">
            <fn:map>
              <fn:string key="name">
                <xsl:value-of select="fn:string[@key='plant_name']"/>
              </fn:string>
            </fn:map>
          </fn:array>
        </xsl:when>
        <xsl:when test="fn:array[@key='plants']/fn:map">
          <fn:array key="organism_sources">
            <xsl:for-each select="fn:array[@key='plants']/fn:map">
              <fn:map>
                <fn:string key="name">
                  <xsl:value-of select="fn:string[@key='name']"/>
                </fn:string>
              </fn:map>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="organism_sources"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">Dr. Duke's</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:string[@key='compound_id']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

  <!-- ============================================================
       Utility Functions
       ============================================================ -->

  <xsl:function name="local:is-null" as="xs:boolean">
    <xsl:param name="value" as="xs:string?"/>
    <xsl:sequence select="not($value) or $value = '' or $value = '-' or $value = 'NA' or $value = 'N/A' or $value = 'null' or $value = '---'"/>
  </xsl:function>

</xsl:stylesheet>
