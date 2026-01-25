<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  Natural Medicines -> Unified Dietary Supplements Schema Mapping
  ============================================================
  Source: ./schema.json
  Target: ../schema.json (Unified Dietary Supplements Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────┬─────────────────┬──────────────┐
  │ Source Field            │ Target Field    │ Transform    │
  ├─────────────────────────┼─────────────────┼──────────────┤
  │ monograph_id            │ ingredient_id   │ Direct       │
  │ name                    │ ingredient_name │ Direct       │
  │ scientific_names        │ scientific_names│ Direct       │
  │ safety_grade            │ safety_rating   │ Direct       │
  │ effectiveness_ratings   │ efficacy        │ Restructure  │
  │ drug_interactions       │ interactions    │ Restructure  │
  │ (computed)              │ _source         │ Metadata     │
  └─────────────────────────┴─────────────────┴──────────────┘

  NOTES:
  - Natural Medicines is monograph-focused (ingredients, not products)
  - Maps to unified ingredient/efficacy schema
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
      <!-- Primary identifier -->
      <fn:string key="ingredient_id">
        <xsl:value-of select="fn:string[@key='monograph_id']"/>
      </fn:string>

      <!-- Ingredient name -->
      <fn:string key="ingredient_name">
        <xsl:value-of select="fn:string[@key='name']"/>
      </fn:string>

      <!-- Scientific names -->
      <xsl:if test="fn:array[@key='scientific_names']">
        <fn:array key="scientific_names">
          <xsl:for-each select="fn:array[@key='scientific_names']/fn:string">
            <fn:string><xsl:value-of select="."/></fn:string>
          </xsl:for-each>
        </fn:array>
      </xsl:if>

      <!-- Common names -->
      <xsl:if test="fn:array[@key='common_names']">
        <fn:array key="common_names">
          <xsl:for-each select="fn:array[@key='common_names']/fn:string">
            <fn:string><xsl:value-of select="."/></fn:string>
          </xsl:for-each>
        </fn:array>
      </xsl:if>

      <!-- Safety rating -->
      <fn:string key="safety_rating">
        <xsl:value-of select="fn:string[@key='safety_grade']"/>
      </fn:string>

      <!-- Efficacy ratings by condition -->
      <xsl:if test="fn:array[@key='effectiveness_ratings']">
        <fn:array key="efficacy">
          <xsl:for-each select="fn:array[@key='effectiveness_ratings']/fn:map">
            <fn:map>
              <fn:string key="condition">
                <xsl:value-of select="fn:string[@key='condition']"/>
              </fn:string>
              <fn:string key="rating">
                <xsl:value-of select="fn:string[@key='rating']"/>
              </fn:string>
              <fn:string key="evidence_level">
                <xsl:value-of select="fn:string[@key='evidence_level']"/>
              </fn:string>
              <xsl:if test="fn:string[@key='typical_dose']">
                <fn:string key="dose">
                  <xsl:value-of select="fn:string[@key='typical_dose']"/>
                </fn:string>
              </xsl:if>
            </fn:map>
          </xsl:for-each>
        </fn:array>
      </xsl:if>

      <!-- Drug interactions -->
      <xsl:if test="fn:array[@key='drug_interactions']">
        <fn:array key="interactions">
          <xsl:for-each select="fn:array[@key='drug_interactions']/fn:map">
            <fn:map>
              <fn:string key="drug">
                <xsl:value-of select="fn:string[@key='drug']"/>
              </fn:string>
              <fn:string key="severity">
                <xsl:value-of select="fn:string[@key='severity']"/>
              </fn:string>
              <xsl:if test="fn:string[@key='mechanism']">
                <fn:string key="mechanism">
                  <xsl:value-of select="fn:string[@key='mechanism']"/>
                </fn:string>
              </xsl:if>
            </fn:map>
          </xsl:for-each>
        </fn:array>
      </xsl:if>

      <!-- Source metadata -->
      <fn:map key="_source">
        <fn:string key="database">Natural Medicines</fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:string[@key='monograph_id']"/>
        </fn:string>
        <fn:string key="record_type">monograph</fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

</xsl:stylesheet>
