<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  Exposome-Explorer -> Unified Metabolomics Schema Mapping
  ============================================================
  Source: ./schema.json
  Target: ../schema.json (Unified Metabolomics Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────┬─────────────────┬──────────────┐
  │ Source Field            │ Target Field    │ Transform    │
  ├─────────────────────────┼─────────────────┼──────────────┤
  │ biomarker.biomarker_name│ metabolite_name │ Direct       │
  │ biomarker.hmdb_id       │ hmdb_id         │ Direct       │
  │ exposure                │ association     │ Restructure  │
  │ specimen_type           │ biofluid        │ Direct       │
  │ effect_size             │ correlation     │ Direct       │
  │ (computed)              │ _source         │ Metadata     │
  └─────────────────────────┴─────────────────┴──────────────┘

  NOTES:
  - Exposome-Explorer focuses on biomarker-exposure associations
  - Maps biomarkers to unified metabolite schema
-->
<xsl:stylesheet version="3.0"
    xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
    xmlns:xs="http://www.w3.org/2001/XMLSchema"
    xmlns:fn="http://www.w3.org/2005/xpath-functions"
    xmlns:map="http://www.w3.org/2005/xpath-functions/map"
    xmlns:array="http://www.w3.org/2005/xpath-functions/array"
    exclude-result-prefixes="xs fn map array">

  <xsl:output method="json" indent="yes"/>

  <xsl:template name="xsl:initial-template">
    <xsl:param name="input" as="xs:string"/>
    <xsl:variable name="source" select="json-to-xml($input)"/>
    <xsl:apply-templates select="$source/fn:map"/>
  </xsl:template>

  <xsl:template match="fn:map">
    <fn:map>
      <!-- Record ID -->
      <fn:number key="record_id">
        <xsl:value-of select="fn:number[@key='association_id']"/>
      </fn:number>

      <!-- Metabolite/biomarker information -->
      <fn:map key="metabolite">
        <fn:string key="name">
          <xsl:value-of select="fn:map[@key='biomarker']/fn:string[@key='biomarker_name']"/>
        </fn:string>
        <xsl:if test="fn:map[@key='biomarker']/fn:string[@key='hmdb_id']">
          <fn:string key="hmdb_id">
            <xsl:value-of select="fn:map[@key='biomarker']/fn:string[@key='hmdb_id']"/>
          </fn:string>
        </xsl:if>
        <xsl:if test="fn:map[@key='biomarker']/fn:number[@key='pubchem_cid']">
          <fn:number key="pubchem_cid">
            <xsl:value-of select="fn:map[@key='biomarker']/fn:number[@key='pubchem_cid']"/>
          </fn:number>
        </xsl:if>
        <xsl:if test="fn:map[@key='biomarker']/fn:string[@key='biomarker_type']">
          <fn:string key="type">
            <xsl:value-of select="fn:map[@key='biomarker']/fn:string[@key='biomarker_type']"/>
          </fn:string>
        </xsl:if>
      </fn:map>

      <!-- Specimen/biofluid -->
      <fn:string key="biofluid">
        <xsl:value-of select="fn:string[@key='specimen_type']"/>
      </fn:string>

      <!-- Exposure association -->
      <fn:map key="exposure">
        <fn:string key="name">
          <xsl:value-of select="fn:map[@key='exposure']/fn:string[@key='exposure_name']"/>
        </fn:string>
        <fn:string key="category">
          <xsl:value-of select="fn:map[@key='exposure']/fn:string[@key='exposure_category']"/>
        </fn:string>
      </fn:map>

      <!-- Association statistics -->
      <fn:map key="association">
        <fn:string key="type">
          <xsl:value-of select="fn:string[@key='association_type']"/>
        </fn:string>
        <xsl:if test="fn:number[@key='effect_size']">
          <fn:number key="effect_size">
            <xsl:value-of select="fn:number[@key='effect_size']"/>
          </fn:number>
        </xsl:if>
        <xsl:if test="fn:string[@key='effect_unit']">
          <fn:string key="effect_unit">
            <xsl:value-of select="fn:string[@key='effect_unit']"/>
          </fn:string>
        </xsl:if>
        <xsl:if test="fn:number[@key='p_value']">
          <fn:number key="p_value">
            <xsl:value-of select="fn:number[@key='p_value']"/>
          </fn:number>
        </xsl:if>
        <xsl:if test="fn:number[@key='sample_size']">
          <fn:number key="sample_size">
            <xsl:value-of select="fn:number[@key='sample_size']"/>
          </fn:number>
        </xsl:if>
      </fn:map>

      <!-- Reference -->
      <xsl:if test="fn:map[@key='reference']/fn:number[@key='pubmed_id']">
        <fn:number key="pubmed_id">
          <xsl:value-of select="fn:map[@key='reference']/fn:number[@key='pubmed_id']"/>
        </fn:number>
      </xsl:if>

      <!-- Source metadata -->
      <fn:map key="_source">
        <fn:string key="database">Exposome-Explorer</fn:string>
        <fn:number key="original_id">
          <xsl:value-of select="fn:number[@key='association_id']"/>
        </fn:number>
      </fn:map>

    </fn:map>
  </xsl:template>

</xsl:stylesheet>
