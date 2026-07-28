---
layout: archive
title: "CMS Notes"
permalink: /notes/
author_profile: true
---

{% include base_path %}

{% for note in site.notes %}
  {% include archive-single.html %}
{% endfor %}