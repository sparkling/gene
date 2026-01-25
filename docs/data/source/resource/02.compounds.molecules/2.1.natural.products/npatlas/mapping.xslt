<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  NPAtlas -> Natural Products Unified Schema Mapping
  ============================================================
  Source: ./schema.md (Natural Products Atlas)
  Target: ../schema.json (2.1 Natural Products Unified Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────┬──────────────┐
  │ Source Field                │ Target Field        │ Transform    │
  ├─────────────────────────────┼─────────────────────┼──────────────┤
  │ npa_id                      │ npa_id              │ Direct       │
  │ npa_id                      │ compound_id         │ Direct       │
  │ compound_name               │ name                │ Direct       │
  │ smiles                      │ canonical_smiles    │ Direct       │
  │ inchi                       │ inchi               │ Direct       │
  │ inchi_key                   │ inchi_key           │ Direct       │
  │ molecular_formula           │ molecular_formula   │ Direct       │
  │ molecular_weight            │ molecular_weight    │ Direct       │
  │ source_organism             │ organism_sources    │ Object array │
  │ compound_class              │ compound_classes    │ Array        │
  │ origin_type                 │ origin_type         │ Enum         │
  │ cluster_type                │ cluster_type        │ Direct       │
  │ isolation_source            │ isolation_source    │ Direct       │
  │ references                  │ references          │ Object array │
  └─────────────────────────────┴─────────────────────┴──────────────┘

  NULL HANDLING:
  - All string fields: "" or "-" → null

  NOTES:
  - NPAtlas uses NPA prefix for compound IDs (e.g., NPA012345)
  - Focuses on microbial natural products (bacterial, fungal, marine)
  - Includes biosynthetic gene cluster information
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
      <!-- ========== IDENTIFIERS ========== -->

      <!-- NPAtlas compound ID -->
      <fn:string key="npa_id">
        <xsl:value-of select="fn:string[@key='npa_id']"/>
      </fn:string>

      <fn:string key="compound_id">
        <xsl:value-of select="fn:string[@key='npa_id']"/>
      </fn:string>

      <!-- ========== BASIC PROPERTIES ========== -->

      <!-- Compound name -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='compound_name'])">
          <fn:null key="name"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="name">
            <xsl:value-of select="fn:string[@key='compound_name']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== CHEMICAL STRUCTURE ========== -->

      <!-- Canonical SMILES -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='smiles'])">
          <fn:null key="canonical_smiles"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="canonical_smiles">
            <xsl:value-of select="fn:string[@key='smiles']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- InChI -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='inchi'])">
          <fn:null key="inchi"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="inchi">
            <xsl:value-of select="fn:string[@key='inchi']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- InChI Key -->
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

      <!-- ========== ORIGIN/SOURCE CLASSIFICATION ========== -->

      <!-- Origin type (bacterial, fungal, marine) -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='origin_type'])">
          <fn:null key="origin_type"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="origin_type">
            <xsl:value-of select="fn:string[@key='origin_type']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Biosynthetic gene cluster type -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='cluster_type'])">
          <fn:null key="cluster_type"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="cluster_type">
            <xsl:value-of select="fn:string[@key='cluster_type']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Isolation source (marine, soil, endophytic, etc.) -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='isolation_source'])">
          <fn:null key="isolation_source"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="isolation_source">
            <xsl:value-of select="fn:string[@key='isolation_source']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== COMPOUND CLASSIFICATION ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='compound_class']/fn:string">
          <fn:array key="compound_classes">
            <xsl:for-each select="fn:array[@key='compound_class']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="compound_classes"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== ORGANISM SOURCES ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='source_organisms']/fn:map">
          <fn:array key="organism_sources">
            <xsl:for-each select="fn:array[@key='source_organisms']/fn:map">
              <fn:map>
                <fn:string key="name">
                  <xsl:value-of select="fn:string[@key='name']"/>
                </fn:string>
                <xsl:if test="fn:number[@key='ncbi_taxon_id']">
                  <fn:number key="ncbi_taxon_id">
                    <xsl:value-of select="fn:number[@key='ncbi_taxon_id']"/>
                  </fn:number>
                </xsl:if>
              </fn:map>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:when test="fn:string[@key='source_organism'] and not(local:is-null(fn:string[@key='source_organism']))">
          <fn:array key="organism_sources">
            <fn:map>
              <fn:string key="name">
                <xsl:value-of select="fn:string[@key='source_organism']"/>
              </fn:string>
            </fn:map>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="organism_sources"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== REFERENCES ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='references']/fn:map">
          <fn:array key="references">
            <xsl:for-each select="fn:array[@key='references']/fn:map">
              <fn:map>
                <xsl:if test="fn:string[@key='doi']">
                  <fn:string key="doi">
                    <xsl:value-of select="fn:string[@key='doi']"/>
                  </fn:string>
                </xsl:if>
                <xsl:if test="fn:number[@key='pmid']">
                  <fn:number key="pmid">
                    <xsl:value-of select="fn:number[@key='pmid']"/>
                  </fn:number>
                </xsl:if>
              </fn:map>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="references"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">NPAtlas</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:string[@key='npa_id']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

  <!-- ============================================================
       Utility Functions
       ============================================================ -->

  <xsl:function name="local:is-null" as="xs:boolean">
    <xsl:param name="value" as="xs:string?"/>
    <xsl:sequence select="not($value) or $value = '' or $value = '-' or $value = 'NA' or $value = 'N/A' or $value = 'null'"/>
  </xsl:function>

</xsl:stylesheet>
