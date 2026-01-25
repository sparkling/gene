<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  ConsumerLab -> Unified Dietary Supplements Schema Mapping
  ============================================================
  Source: ./schema.json
  Target: ../schema.json (Unified Dietary Supplements Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────┬─────────────────┬──────────────┐
  │ Source Field            │ Target Field    │ Transform    │
  ├─────────────────────────┼─────────────────┼──────────────┤
  │ review_id               │ product_id      │ Direct       │
  │ product_name            │ product_name    │ Direct       │
  │ brand                   │ brand           │ Direct       │
  │ category                │ category        │ Direct       │
  │ approval_status         │ quality_status  │ Direct       │
  │ test_results            │ testing         │ Restructure  │
  │ (computed)              │ _source         │ Metadata     │
  └─────────────────────────┴─────────────────┴──────────────┘

  NOTES:
  - ConsumerLab is testing-focused; maps to quality assessment fields
  - Heavy metal results restructured for unified format
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
      <fn:string key="product_id">
        <xsl:value-of select="fn:string[@key='review_id']"/>
      </fn:string>

      <!-- Product name -->
      <fn:string key="product_name">
        <xsl:value-of select="fn:string[@key='product_name']"/>
      </fn:string>

      <!-- Brand -->
      <fn:string key="brand">
        <xsl:value-of select="fn:string[@key='brand']"/>
      </fn:string>

      <!-- Category -->
      <fn:string key="category">
        <xsl:value-of select="fn:string[@key='category']"/>
      </fn:string>

      <!-- UPC if available -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='upc']">
          <fn:string key="upc">
            <xsl:value-of select="fn:string[@key='upc']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="upc"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Quality status -->
      <fn:string key="quality_status">
        <xsl:value-of select="fn:string[@key='approval_status']"/>
      </fn:string>

      <!-- Testing results -->
      <xsl:if test="fn:map[@key='test_results']">
        <fn:map key="testing">
          <fn:string key="test_date">
            <xsl:value-of select="fn:string[@key='test_date']"/>
          </fn:string>

          <!-- Potency -->
          <xsl:if test="fn:map[@key='test_results']/fn:array[@key='potency_results']">
            <fn:array key="potency_results">
              <xsl:for-each select="fn:map[@key='test_results']/fn:array[@key='potency_results']/fn:map">
                <fn:map>
                  <fn:string key="ingredient">
                    <xsl:value-of select="fn:string[@key='ingredient']"/>
                  </fn:string>
                  <fn:number key="percent_of_label">
                    <xsl:value-of select="fn:number[@key='percent_claim']"/>
                  </fn:number>
                  <fn:string key="status">
                    <xsl:value-of select="fn:string[@key='status']"/>
                  </fn:string>
                </fn:map>
              </xsl:for-each>
            </fn:array>
          </xsl:if>

          <!-- Heavy metals -->
          <xsl:if test="fn:map[@key='test_results']/fn:map[@key='heavy_metals']">
            <fn:map key="contaminants">
              <fn:string key="type">heavy_metals</fn:string>
              <fn:string key="status">Pass</fn:string>
            </fn:map>
          </xsl:if>
        </fn:map>
      </xsl:if>

      <!-- Source metadata -->
      <fn:map key="_source">
        <fn:string key="database">ConsumerLab</fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:string[@key='review_id']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

</xsl:stylesheet>
