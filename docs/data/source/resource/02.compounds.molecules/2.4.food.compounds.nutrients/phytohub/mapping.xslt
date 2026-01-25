<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  PhytoHub -> Food Compounds Unified Schema Mapping
  ============================================================
  Source: ./schema.md (PhytoHub Dietary Phytochemicals Database)
  Target: ../schema.json (2.4 Food Compounds and Nutrients Unified Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────┬──────────────┐
  │ Source Field                │ Target Field        │ Transform    │
  ├─────────────────────────────┼─────────────────────┼──────────────┤
  │ phytohub_id                 │ phytohub_id         │ Direct       │
  │ phytohub_id                 │ compound_id         │ Direct       │
  │ name                        │ name                │ Direct       │
  │ smiles                      │ smiles              │ Direct       │
  │ inchi_key                   │ inchi_key           │ Direct       │
  │ molecular_weight            │ molecular_weight    │ Direct       │
  │ compound_type               │ compound_type       │ Direct       │
  │ parent_id                   │ parent_id           │ Direct       │
  │ transformation              │ transformation      │ Direct       │
  │ ms_spectra                  │ ms_spectra          │ Object array │
  │ external_ids                │ external_ids        │ Object array │
  │ food_content                │ food_sources        │ Object array │
  └─────────────────────────────┴─────────────────────┴──────────────┘

  NULL HANDLING:
  - All string fields: "" → null

  NOTES:
  - PhytoHub uses PHUB prefix for IDs (e.g., PHUB000001)
  - Includes parent-metabolite relationships
  - MS/MS reference spectra for metabolomics
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

  <xsl:template name="xsl:initial-template">
    <xsl:param name="input" as="xs:string"/>
    <xsl:variable name="source" select="json-to-xml($input)"/>
    <xsl:apply-templates select="$source/fn:map"/>
  </xsl:template>

  <xsl:template match="fn:map">
    <fn:map>
      <!-- ========== IDENTIFIERS ========== -->

      <fn:string key="phytohub_id">
        <xsl:value-of select="fn:string[@key='phytohub_id']"/>
      </fn:string>

      <fn:string key="compound_id">
        <xsl:value-of select="fn:string[@key='phytohub_id']"/>
      </fn:string>

      <!-- ========== BASIC PROPERTIES ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='name'])">
          <fn:null key="name"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="name">
            <xsl:value-of select="fn:string[@key='name']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== CHEMICAL STRUCTURE ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='smiles'])">
          <fn:null key="smiles"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="smiles">
            <xsl:value-of select="fn:string[@key='smiles']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='inchi_key'])">
          <fn:null key="inchi_key"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="inchi_key">
            <xsl:value-of select="fn:string[@key='inchi_key']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

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

      <!-- ========== COMPOUND TYPE / METABOLITE RELATIONSHIP ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='compound_type'])">
          <fn:null key="compound_type"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="compound_type">
            <xsl:value-of select="fn:string[@key='compound_type']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='parent_id'])">
          <fn:null key="parent_id"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="parent_id">
            <xsl:value-of select="fn:string[@key='parent_id']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='transformation'])">
          <fn:null key="transformation"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="transformation">
            <xsl:value-of select="fn:string[@key='transformation']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== MS/MS SPECTRA ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='ms_spectra']/fn:map">
          <fn:array key="ms_spectra">
            <xsl:for-each select="fn:array[@key='ms_spectra']/fn:map">
              <fn:map>
                <fn:number key="mz">
                  <xsl:value-of select="fn:number[@key='mz']"/>
                </fn:number>
                <fn:number key="intensity">
                  <xsl:value-of select="fn:number[@key='intensity']"/>
                </fn:number>
              </fn:map>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="ms_spectra"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== EXTERNAL IDS ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='external_ids']/fn:map">
          <fn:array key="external_ids">
            <xsl:for-each select="fn:array[@key='external_ids']/fn:map">
              <fn:map>
                <fn:string key="database">
                  <xsl:value-of select="fn:string[@key='database']"/>
                </fn:string>
                <fn:string key="id">
                  <xsl:value-of select="fn:string[@key='id']"/>
                </fn:string>
              </fn:map>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="external_ids"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== FOOD SOURCES ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='food_content']/fn:map">
          <fn:array key="food_sources">
            <xsl:for-each select="fn:array[@key='food_content']/fn:map">
              <fn:map>
                <fn:string key="food_name">
                  <xsl:value-of select="fn:string[@key='food_name']"/>
                </fn:string>
                <xsl:if test="fn:number[@key='content_value']">
                  <fn:number key="content_value">
                    <xsl:value-of select="fn:number[@key='content_value']"/>
                  </fn:number>
                </xsl:if>
                <xsl:if test="fn:string[@key='content_unit']">
                  <fn:string key="content_unit">
                    <xsl:value-of select="fn:string[@key='content_unit']"/>
                  </fn:string>
                </xsl:if>
              </fn:map>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="food_sources"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">PhytoHub</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:string[@key='phytohub_id']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

  <xsl:function name="local:is-null" as="xs:boolean">
    <xsl:param name="value" as="xs:string?"/>
    <xsl:sequence select="not($value) or $value = '' or $value = '-' or $value = 'NA' or $value = 'N/A' or $value = 'null'"/>
  </xsl:function>

</xsl:stylesheet>
